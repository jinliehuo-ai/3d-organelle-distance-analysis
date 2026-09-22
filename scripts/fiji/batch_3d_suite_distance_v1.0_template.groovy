// PUBLIC TEMPLATE | analysis baseline 1.0 | 2026-09-22
// Copy this file before changing an algorithm, threshold, inclusion rule, or output field.
import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.measure.Calibration
import ij.plugin.GaussianBlur3D
import ij.process.ByteProcessor
import ij.process.ShortProcessor
import ij.process.ColorProcessor
import ij.process.ImageProcessor
import ij.process.StackStatistics
import ij.process.AutoThresholder
import ij.io.FileSaver
import loci.plugins.BF
import mcib3d.image3d.ImageHandler
import mcib3d.image3d.ImageInt
import mcib3d.image3d.ImageLabeller
import mcib3d.geom.Object3D
import mcib3d.geom.Objects3DPopulation
import mcib3d.geom.Voxel3D

// Batch 3D Suite analysis: one nucleus and one spot per single-cell image.
//
// The 3D Manager GUI is deliberately NOT used. Its "Distances" button reads
// javax.swing.JList.getSelectedIndices() and indexes the object table with the
// result, so programmatically added objects and a reused manager window make
// the indices drift out of range (ArrayIndexOutOfBoundsException) and leave
// stale objects from the previous image in the table. Every column of the
// Manager's "3D Distance" window is reproduced below through the mcib3d API,
// which is deterministic, headless-safe and always emits exactly N*(N-1) rows.
// Public template paths are relative to -DprojectDir or the current directory.
// Example: fiji --headless --run '...groovy', 'projectDir=/path/to/project'
def projectDir = new File(System.getProperty("projectDir", "."))
def inputDir = new File(projectDir, "data/raw")
def outputDir = new File(projectDir, "results")
outputDir.mkdirs()
new File(outputDir, "all_files_summary.csv").delete()
new File(outputDir, "all_files_summary.xls").delete()
// Per-image output is filed under the cell type, the leading letters of the name
// (SampleA1-2 -> SampleA), so one run of many groups does not leave a flat directory
// that has to be sorted by eye. Outputs that span groups stay at the top.
inputDir.eachFileRecurse { file ->
    if (!file.isFile()) {
        return
    }
    String oldStem = file.name.replaceFirst(/(?i)\.(nd2|tif|tiff)$/, "")
    def staleDirs = [outputDir]
    outputDir.eachFile { entry -> if (entry.isDirectory()) staleDirs << entry }
    staleDirs.each { dir ->
        new File(dir, oldStem + "_spot_to_nucleus_distance.csv").delete()
        new File(dir, oldStem + "_manager_distances.csv").delete()
        new File(dir, oldStem + "_spot_to_nucleus_distance.xls").delete()
        new File(dir, oldStem + "_manager_distances.xls").delete()
        new File(dir, oldStem + "_ERROR.txt").delete()
        new File(dir, oldStem + "_metadata.txt").delete()
        new File(dir, oldStem + "_nuclei_labels.tif").delete()
        new File(dir, oldStem + "_spots_labels.tif").delete()
    }
    // Overlays wherever an earlier grouping rule put them.
    staleDirs.each { dir -> new File(new File(dir, "qc_overlays"), "qc_" + oldStem + ".png").delete() }
}

final int NUCLEUS_CHANNEL = 1
final int SPOT_CHANNEL = 3
final double NUCLEUS_SIGMA_XY = 1.0
final double NUCLEUS_SIGMA_Z = 1.0
final double SPOT_SIGMA_XY = 0.5
final double SPOT_SIGMA_Z = 0.5
final int MIN_NUCLEUS_VOXELS = 10
final int MIN_SPOT_VOXELS = 20
// The hard reject for an object that is not a centrosome at all. It is NOT a
// size prior on the centrosome: MaxEntropy sits low relative to a punctum that
// saturates the detector, so the segmented blob is the centrosome plus the halo
// the threshold could not separate from it, and a real one reaches 6684 voxels
// in two crops. At a 5000-voxel cap those two were rejected outright and the
// fallback was a 134-voxel sliver against the image edge, 4.58 um from the
// nucleus, instead of the punctum at 2.80 um. Raising the threshold on such a
// blob shrinks it smoothly to a compact core without ever splitting it and moves
// its centroid by less than a voxel, which is what says it is one object with a
// halo rather than an aggregate. So this only has to stop the whole-cytoplasm
// blob an over-permissive threshold produces, which is an order of magnitude
// larger again. PLAUSIBLE_SPOT_VOXELS stays the soft warning about size.
final int MAX_SPOT_VOXELS = 10000
// Otsu and Default (IsoData) assume a bimodal foreground/background histogram,
// which holds for DAPI but not for a sparse punctate marker: on C3 they cut at
// ~230/4095 and merge the diffuse cytoplasmic signal into one 150k-voxel object
// that overlaps the nucleus, collapsing the measured distance to ~0.5 um.
// MaxEntropy cuts at ~1180 and isolates the ~300-voxel centrosome, reproducing
// the manual measurement (one check crop: 2.53 um vs 2.0-2.5 um measured by hand).
final AutoThresholder.Method NUCLEUS_THRESHOLD_METHOD = AutoThresholder.Method.Otsu
final AutoThresholder.Method SPOT_THRESHOLD_METHOD = AutoThresholder.Method.MaxEntropy
// "intensity" keeps the brightest candidate, "volume" keeps the largest one.
// A gamma-tubulin punctum is defined by brightness, not by size.
final String SPOT_SELECT_BY = "intensity"
// Axial resolution of the acquisition, FWHM_z ~ 1.4 * lambda * n / NA^2; for
// 520 nm through a 1.4 NA oil objective that is ~0.56 um. A spot whose peak
// lies closer to the nuclear surface than this is docked on the envelope: the
// separation, and even its sign, is below what the microscope can resolve, so
// reporting a signed sub-resolution number would be false precision. Set this
// from the objective actually used.
final double AXIAL_RESOLUTION_UM = 0.5
// How much of a small sphere centred on the spot is nucleus. This is the
// geometry the distance is trying to describe, read directly and without
// depending on a ray direction: ~0 means the spot is free in cytoplasm, ~0.5
// means it is docked on the envelope with the nucleus filling one side, and
// ~1 means the mask has closed around it. DAPI intensity at the spot cannot
// substitute for this, because the axial PSF spills nuclear signal about
// 0.3 um past the envelope, so a correctly docked spot also reads bright.
final double SHELL_RADIUS_UM = 1.0
final double DOCKED_SHELL_FRACTION = 0.30
final double ENGULFED_SHELL_FRACTION = 0.70
// The same fraction recomputed against a stricter threshold. A spot sitting on
// a steep nuclear edge barely moves between the two; one sitting where DAPI
// hovers at the cut swings, and its distance is a threshold artefact.
final double THRESHOLD_SENSITIVE_DELTA = 0.10
// A 12-bit detector saturates at 4095, so two genuine puncta in one crop can
// land on exactly the same peak value and the ranking then turns on whatever
// tiebreak is applied, which carries no biological meaning. Only candidates
// this close in brightness count as tied, and they are settled on proximity
// instead: a cell's own centrosome sits beside its own nucleus, while an
// equally bright punctum several microns away usually belongs to the
// neighbouring cell that a single-cell crop clipped into frame. The cut is
// deliberately near 1: a candidate that is merely dimmer, not tied, must still
// lose on brightness, or a small dim blob hugging the nucleus would win.
final double SPOT_BRIGHTNESS_TOLERANCE = 0.95
// Beyond this a punctum is too far from the nucleus to be that cell's own
// centrosome; the row is kept but flagged rather than quietly averaged in.
final double IMPLAUSIBLE_DISTANCE_UM = 5.0
// A centrosome is roughly 0.5-1 um across; at this voxel size a 2 um sphere is
// already ~430 voxels. Past this the object is a merged pair, an aggregate or a
// stretch of cytoplasmic signal rather than one centrosome, so the row is
// flagged. MAX_SPOT_VOXELS stays the hard reject; this is the softer warning.
final int PLAUSIBLE_SPOT_VOXELS = 1000
// A duplicated centrosome pair sits well within a micron or two of itself, so a
// second candidate this close is the partner rather than a different structure.
final double PAIR_SEPARATION_UM = 2.5
// A nucleus that runs off the first or last slice, or out of the XY field, is
// cut by a flat plane that is an artefact of the acquisition, not its envelope.
// mcib3d builds an object's contour by testing 6-neighbours with everything
// outside the image counted as background, so those cut faces enter the contour
// and distCenterBorderUnit is free to return the distance to the crop instead of
// to the nuclear surface. Surface voxels that are surface only because the stack
// ended are dropped before the nearest-surface search below.
final double SURFACE_CLIP_TOLERANCE_UM = 1e-9
// How much of the nuclear surface may be that artificial plane before the whole
// measurement is called into question. Below this the envelope is still mostly
// present and the nearest-surface search has real surface to find; above it the
// object in the stack is a slab of a nucleus and even an uncontaminated distance
// describes a shape that was only partly acquired.
final double TRUNCATED_SURFACE_FRACTION = 0.15
// A shell sample landing outside the acquired volume counts against the
// fraction: nothing was recorded there, so it is read as "not nucleus". An
// earlier version dropped those directions instead and renormalised over the
// imaged part of the sphere, which divided by a shrinking denominator exactly
// where the missing part was the cytoplasmic side and read two docked spots as
// engulfed. Past this much of the shell missing, the row is flagged and
// shell_upper says how far above the reported fraction the truth could lie.
final double SHELL_OUT_OF_VOLUME_FRACTION = 0.10

File summaryFile = new File(outputDir, "all_files_summary.xls")
summaryFile.text = "file_name\tdistance_to_nucleus_surface_um\tnegative_distance\tdistance_from_peak_um\tsurface_to_surface_um\tqc_status\tshell_in_nucleus\tshell_strict\tdapi_at_spot_ratio\tspot_voxels\tspot_n\tnuclei_n\tpair_separation_um\tsecond_spot_um\tspot_in_nucleus_pct\tspot_z_clipped\tnucleus_z_clipped\tnucleus_surface_clipped_frac\tshell_out_of_volume\tshell_lower\tshell_upper\tspot_id\tnucleus_id\tspot_x_um\tspot_y_um\tspot_z_um\n"
File logFile = new File(outputDir, "batch_log.txt")
logFile.text = "Input folder: ${inputDir.absolutePath}\n"

// The filesystem hands files back in arbitrary order, which scatters the summary
// rows. Sorting here puts both the summary and the log in name order.
def imageFiles = []
inputDir.eachFileRecurse { candidate ->
    if (candidate.isFile() && (candidate.name.toLowerCase().endsWith(".nd2") || candidate.name.toLowerCase().endsWith(".tif") || candidate.name.toLowerCase().endsWith(".tiff"))) {
        imageFiles << candidate
    }
}
// sort(false) returns a new list rather than relying on in-place mutation,
// which is why the result has to be assigned back.
imageFiles = imageFiles.sort(false) { naturalKey(it.name) }

imageFiles.each { file ->

    String stem = file.name.replaceFirst(/(?i)\.(nd2|tif|tiff)$/, "")
    File groupDir = new File(outputDir, groupOf(stem))
    groupDir.mkdirs()
    // Every run leaves a rendered check of what was actually segmented, so a
    // wrong object never has to be caught by reading numbers alone.
    File overlayDir = new File(groupDir, "qc_overlays")
    overlayDir.mkdirs()
    IJ.log("Processing: " + file.name)
    logFile << "Processing: ${file.name}\n"
    try {
    ImagePlus imp = file.name.toLowerCase().endsWith(".nd2") ? BF.openImagePlus(file.absolutePath)[0] : IJ.openImage(file.absolutePath)
    if (imp == null) {
        IJ.log("Skipped: could not open " + file.name)
        return
    }
    if (imp.getNChannels() < SPOT_CHANNEL) {
        IJ.log("Skipped: fewer than 3 channels in " + file.name)
        imp.close()
        return
    }

    Calibration cal = imp.getCalibration().copy()
    double sx = cal.pixelWidth
    double sy = cal.pixelHeight
    double sz = cal.pixelDepth
    if (!(sx > 0) || !(sy > 0) || !(sz > 0)) {
        throw new IllegalArgumentException("Missing voxel size in " + file.name)
    }
    if (Math.abs(sx - sy) > 1e-9) {
        throw new IllegalArgumentException("Anisotropic XY voxels in ${file.name}: mcib3d takes a single XY scale")
    }

    ImagePlus nucleiRaw = extractChannel(imp, NUCLEUS_CHANNEL)
    ImagePlus spotsRaw = extractChannel(imp, SPOT_CHANNEL)
    // Intensities are read back from the unfiltered spot channel so that the
    // ranking is not biased by the blur.
    ImagePlus spotDisplay = extractChannel(imp, SPOT_CHANNEL)
    ImagePlus nucleusDisplay = extractChannel(imp, NUCLEUS_CHANNEL)
    ImageHandler spotIntensity = ImageHandler.wrap(spotDisplay)
    ImageHandler nucleusIntensity = ImageHandler.wrap(nucleusDisplay)
    GaussianBlur3D.blur(nucleiRaw, NUCLEUS_SIGMA_XY, NUCLEUS_SIGMA_XY, NUCLEUS_SIGMA_Z)
    GaussianBlur3D.blur(spotsRaw, SPOT_SIGMA_XY, SPOT_SIGMA_XY, SPOT_SIGMA_Z)

    ImagePlus nucleiMask = thresholdStack(nucleiRaw, NUCLEUS_THRESHOLD_METHOD, 0.0)
    StackStatistics spotStatistics = new StackStatistics(spotsRaw)
    int spotThresholdBin = new AutoThresholder().getThreshold(SPOT_THRESHOLD_METHOD, intHistogram(spotStatistics.histogram))
    double spotThreshold = spotStatistics.histMin + (spotThresholdBin + 0.5) * spotStatistics.binSize
    ImagePlus spotsMask = thresholdStack(spotsRaw, null, spotThreshold)

    ImageHandler nucleiHandler = ImageHandler.wrap(nucleiMask)
    ImageHandler spotsHandler = ImageHandler.wrap(spotsMask)
    ImageLabeller nucleiLabeller = new ImageLabeller()
    ImageLabeller spotsLabeller = new ImageLabeller()
    nucleiLabeller.setMinSize(MIN_NUCLEUS_VOXELS)
    spotsLabeller.setMinSize(MIN_SPOT_VOXELS)
    ImageInt nucleiLabels = nucleiLabeller.getLabels(nucleiHandler, false)
    ImageInt spotsLabels = spotsLabeller.getLabels(spotsHandler, false)
    Objects3DPopulation nuclei = new Objects3DPopulation(nucleiLabels, 0)
    Objects3DPopulation spots = new Objects3DPopulation(spotsLabels, 0)
    nuclei.setCalibration(sx, sz, cal.unit)
    spots.setCalibration(sx, sz, cal.unit)

    def nucleusCandidates = nuclei.getObjectsList().findAll { it.getVolumePixels() >= MIN_NUCLEUS_VOXELS }
    def spotCandidates = spots.getObjectsList().findAll { it.getVolumePixels() >= MIN_SPOT_VOXELS && it.getVolumePixels() <= MAX_SPOT_VOXELS }
    def oversizedSpots = spots.getObjectsList().findAll { it.getVolumePixels() > MAX_SPOT_VOXELS }
    def nucleiList = nucleusCandidates.isEmpty() ? [] : [nucleusCandidates.max { it.getVolumePixels() }]
    // A 12-bit detector saturates at 4095, so peak intensity alone ties between
    // candidates; mean intensity breaks the tie. A two-argument closure is used
    // as a comparator because a one-argument closure would have to return a
    // single Comparable, and a list of the two keys is not one.
    def brightest = { Object3D a, Object3D b ->
        int byPeak = Double.compare(a.getPixMaxValue(spotIntensity), b.getPixMaxValue(spotIntensity))
        byPeak != 0 ? byPeak : Double.compare(a.getPixMeanValue(spotIntensity), b.getPixMeanValue(spotIntensity))
    }
    // A centrosome is never inside the nucleus, so ranking candidates by
    // brightness alone is not enough: where an intranuclear object in the
    // gamma-tubulin channel outshines the real centrosome, brightness picks the
    // wrong one. The constraint is applied as a filter, and whatever it rejects
    // is written to the log rather than silently dropped.
    def brightestOverall = spotCandidates.isEmpty() ? null : spotCandidates.max(brightest)
    def extranuclear = []
    if (!nucleiList.isEmpty()) {
        Object3D candidateNucleus = nucleiList[0]
        ImagePlus nucleusProbeImage = blankLike(imp)
        ImageHandler nucleusProbe = ImageHandler.wrap(nucleusProbeImage)
        candidateNucleus.draw(nucleusProbe, 255)
        extranuclear = spotCandidates.findAll { Object3D candidate ->
            def candidatePeak = candidate.getPixelMax(spotIntensity)
            shellFractionInMask(nucleusProbe, candidatePeak, SHELL_RADIUS_UM, sx, sz, imp) < ENGULFED_SHELL_FRACTION
        }
        nucleusProbe.closeImagePlus()
    }
    def spotsList = []
    boolean brighterIntranuclearRejected = false
    if (SPOT_SELECT_BY != "intensity") {
        spotsList = spotCandidates.isEmpty() ? [] : [spotCandidates.max { it.getVolumePixels() }]
    } else if (!extranuclear.isEmpty()) {
        Object3D ownNucleus = nucleiList[0]
        double bestPeak = extranuclear.collect { it.getPixMaxValue(spotIntensity) }.max()
        def comparablyBright = extranuclear.findAll {
            it.getPixMaxValue(spotIntensity) >= bestPeak * SPOT_BRIGHTNESS_TOLERANCE
        }
        Object3D chosen = comparablyBright.min { it.distCenterBorderUnit(ownNucleus) }
        if (comparablyBright.size() > 1) {
            logFile << "${file.name}: ${comparablyBright.size()} candidates within ${(int) (SPOT_BRIGHTNESS_TOLERANCE * 100)}% of peak brightness (tied); took the one nearest the nucleus (${String.format('%.3f', chosen.distCenterBorderUnit(ownNucleus))} um) over ${comparablyBright.collect { String.format('%.3f', it.distCenterBorderUnit(ownNucleus)) }.join(', ')} um\n"
        }
        spotsList = [chosen]
        // Compared on intensity rather than object identity: that is the fact
        // worth flagging, and it does not depend on how Object3D defines equality.
        if (brightestOverall != null && brightestOverall.getPixMaxValue(spotIntensity) > chosen.getPixMaxValue(spotIntensity)) {
            brighterIntranuclearRejected = true
            logFile << "${file.name}: rejected a brighter intranuclear candidate (maxI=${brightestOverall.getPixMaxValue(spotIntensity)}, ${brightestOverall.getVolumePixels()} vox) in favour of an extranuclear one (maxI=${chosen.getPixMaxValue(spotIntensity)}, ${chosen.getVolumePixels()} vox)\n"
        }
    } else if (brightestOverall != null) {
        // Every candidate is engulfed; keep the brightest so the row still
        // carries the evidence, and let the QC status say it is unusable.
        spotsList = [brightestOverall]
        logFile << "${file.name}: no extranuclear spot candidate; every candidate sits inside the nucleus mask\n"
    }
    if (!oversizedSpots.isEmpty()) {
        logFile << "${file.name}: WARNING ${oversizedSpots.size()} spot object(s) above ${MAX_SPOT_VOXELS} voxels rejected (largest=${oversizedSpots.max { it.getVolumePixels() }.getVolumePixels()}); the spot threshold may be too permissive\n"
    }
    logFile << "${file.name}: nuclei=${nucleiList.size()}, spots=${spotsList.size()}\n"

    if (nucleiList.isEmpty() || spotsList.isEmpty()) {
        logFile << "${file.name}: skipped, need one nucleus and one spot (nucleus candidates=${nucleusCandidates.size()}, spot candidates=${spotCandidates.size()})\n"
        IJ.log("Skipped: no nucleus/spot pair in " + file.name)
        spotIntensity.closeImagePlus()
        nucleiRaw.close(); spotsRaw.close(); nucleiMask.close(); spotsMask.close(); imp.close()
        return
    }

    Object3D nucleus = nucleiList[0]
    Object3D spot = spotsList[0]
    nucleus.setValue(1)
    spot.setValue(2)
    nucleus.setName("Nucleus_1")
    spot.setName("Spot_1")

    // Object order matches what the 3D Manager would list: nucleus first, then
    // spot, so Obj1/Obj2 in the table keep the meaning they had in the GUI.
    def managed = [nucleus, spot]
    def population = new Objects3DPopulation(managed as Object3D[])
    population.setCalibration(sx, sz, cal.unit)

    // Quality control. distCenterBorderUnit is an UNSIGNED distance to the
    // nearest surface voxel, so a spot buried inside the nucleus mask yields a
    // small positive number that looks like a tight association but is
    // physically meaningless: a centrosome cannot sit inside the nucleus.
    // Overlap is measured directly rather than inferred from bor-bor == 0,
    // which is also produced by two objects merely touching.
    ImagePlus nucleusMaskImage = blankLike(imp)
    ImageHandler nucleusOnly = ImageHandler.wrap(nucleusMaskImage)
    nucleus.draw(nucleusOnly, 255)
    int overlapVoxels = 0
    spot.getVoxels().each { voxel ->
        int vx = (int) voxel.getX(), vy = (int) voxel.getY(), vz = (int) voxel.getZ()
        if (vx >= 0 && vx < imp.getWidth() && vy >= 0 && vy < imp.getHeight() && vz >= 0 && vz < imp.getNSlices()
                && nucleusOnly.getPixel(vx, vy, vz) > 0) {
            overlapVoxels++
        }
    }
    double overlapPercent = 100.0 * overlapVoxels / spot.getVolumePixels()
    // A structure whose voxels reach the first or last slice continues outside
    // the acquired volume, so both its centroid and the nuclear surface facing
    // it are cut off by the crop rather than measured.
    int lastSlice = imp.getNSlices() - 1
    boolean spotZClipped = spot.getZmin() <= 0 || spot.getZmax() >= lastSlice
    boolean nucleusZClipped = nucleus.getZmin() <= 0 || nucleus.getZmax() >= lastSlice

    // The nuclear surface as the distance functions see it, and the same set with
    // the crop faces removed. Every nearest-surface search below is run against
    // both: where they agree, mcib3d's own number is reported unchanged; where
    // the free surface is farther, mcib3d had locked onto the truncation plane
    // and the row is recomputed and flagged.
    ImagePlus spotMaskImage = blankLike(imp)
    ImageHandler spotOnly = ImageHandler.wrap(spotMaskImage)
    spot.draw(spotOnly, 255)
    List<Voxel3D> nucleusSurface = surfaceVoxels(nucleus, nucleusOnly, imp)
    List<Voxel3D> nucleusSurfaceFree = nucleusSurface.findAll { !isTruncationFace(it, nucleusOnly, imp) }
    List<Voxel3D> spotSurface = surfaceVoxels(spot, spotOnly, imp)
    double truncatedSurfaceFraction = nucleusSurface.isEmpty() ? Double.NaN
            : 1.0 - (double) nucleusSurfaceFree.size() / nucleusSurface.size()
    boolean nucleusTruncated = truncatedSurfaceFraction > TRUNCATED_SURFACE_FRACTION
    // A nucleus with no free surface at all would leave nothing to measure to;
    // fall back to the whole contour so the row still carries a number, and let
    // NUCLEUS_TRUNCATED say what it is worth.
    if (nucleusSurfaceFree.isEmpty()) {
        nucleusSurfaceFree = nucleusSurface
    }

    // Placeholder; the status needs the peak distance computed below.
    String qcStatus = "OK"

    if (qcStatus != "OK") {
        logFile << "${file.name}: QC ${qcStatus} (shell in nucleus ${String.format('%.2f', shellFraction)}, strict ${String.format('%.2f', shellFractionStrict)}, DAPI at spot ${String.format('%.2f', dapiRatio)}x nuclear median, peak ${String.format('%+.3f', peakToNucleusBorder)} um from surface, blob ${String.format('%.0f', overlapPercent)}% inside nucleus, spotZ[${spot.getZmin()}..${spot.getZmax()}], nucleusZ[${nucleus.getZmin()}..${nucleus.getZmax()}], lastSlice=${lastSlice})\n"
        IJ.log("QC ${qcStatus} in ${file.name}")
    }

    File managerFile = new File(groupDir, stem + "_manager_distances.xls")
    managerFile.text = "Nb\tObj1\tObj2\tType1\tType2\tLabel1\tLabel2\tcen-cen\tcen-bor\tbor-bor\tradiusCen\texcen\tbor-rad\tperiph\tclosest_cen_i\tclosest_bor_i\tclosest_cen_n\tclosest_bor_n\n"
    distanceRows(managed, population).each { managerFile << it }

    // distCenterBorderUnit measures to the nearest surface voxel from either
    // side, so on its own it reports a spot buried 0.1 um inside the nucleus and
    // one sitting 0.1 um outside it as the same number. The sign is restored
    // from an explicit inside test.
    int centroidX = (int) Math.round(spot.getCenterX())
    int centroidY = (int) Math.round(spot.getCenterY())
    int centroidZ = (int) Math.round(spot.getCenterZ())
    boolean centroidInside = insideMask(nucleusOnly, centroidX, centroidY, centroidZ, imp)
    // The truncation plane is only a problem when it is what the search settles
    // on. Comparing the two searches voxel for voxel decides that per image, so
    // an image whose crop is nowhere near the spot keeps mcib3d's number exactly.
    double centroidToAll = nearestSurfaceUnit(nucleusSurface, spot.getCenterX(), spot.getCenterY(), spot.getCenterZ(), sx, sy, sz)
    double centroidToFree = nearestSurfaceUnit(nucleusSurfaceFree, spot.getCenterX(), spot.getCenterY(), spot.getCenterZ(), sx, sy, sz)
    boolean centroidSurfaceClipped = centroidToFree - centroidToAll > SURFACE_CLIP_TOLERANCE_UM
    double spotCenterToNucleusBorder = (centroidSurfaceClipped ? centroidToFree : spot.distCenterBorderUnit(nucleus)) * (centroidInside ? -1.0 : 1.0)

    // The axial PSF smears a punctum along z, which drags its centroid but not
    // its intensity peak, so the peak is the more stable anchor for the spot.
    Voxel3D peak = brightestPlateauCentre(spot, spotIntensity)
    boolean peakInside = insideMask(nucleusOnly, (int) peak.getX(), (int) peak.getY(), (int) peak.getZ(), imp)
    double peakToAll = nearestSurfaceUnit(nucleusSurface, peak.getX(), peak.getY(), peak.getZ(), sx, sy, sz)
    double peakToFree = nearestSurfaceUnit(nucleusSurfaceFree, peak.getX(), peak.getY(), peak.getZ(), sx, sy, sz)
    boolean peakSurfaceClipped = peakToFree - peakToAll > SURFACE_CLIP_TOLERANCE_UM
    double peakToNucleusBorder = (peakSurfaceClipped ? peakToFree
            : nucleus.distPixelBorderUnit(peak.getX(), peak.getY(), peak.getZ())) * (peakInside ? -1.0 : 1.0)

    double borderToAll = nearestSurfacePair(spotSurface, nucleusSurface, sx, sy, sz)
    double borderToFree = nearestSurfacePair(spotSurface, nucleusSurfaceFree, sx, sy, sz)
    boolean borderSurfaceClipped = borderToFree - borderToAll > SURFACE_CLIP_TOLERANCE_UM
    double surfaceToSurface = borderSurfaceClipped ? borderToFree : spot.distBorderUnit(nucleus)
    boolean surfaceClipped = centroidSurfaceClipped || peakSurfaceClipped || borderSurfaceClipped

    // How much nuclear signal is present where the spot is. Voxels within 3 um
    // of the peak are excluded from the reference so that a spot sitting in an
    // invagination cannot drag down the very median it is compared against.
    List<Double> nucleusReference = []
    nucleus.getVoxels().each { voxel ->
        double dx = (voxel.getX() - peak.getX()) * sx
        double dy = (voxel.getY() - peak.getY()) * sy
        double dz = (voxel.getZ() - peak.getZ()) * sz
        if (Math.sqrt(dx * dx + dy * dy + dz * dz) > 3.0) {
            nucleusReference << (double) nucleusIntensity.getPixel((int) voxel.getX(), (int) voxel.getY(), (int) voxel.getZ())
        }
    }
    nucleusReference.sort()
    double nucleusMedian = nucleusReference.isEmpty() ? Double.NaN : nucleusReference[(int) (nucleusReference.size() / 2)]
    double dapiAtSpot = nucleusIntensity.getPixel((int) peak.getX(), (int) peak.getY(), (int) peak.getZ())
    double dapiRatio = nucleusMedian > 0 ? dapiAtSpot / nucleusMedian : Double.NaN

    // The axially smeared segmented blob can overlap the nucleus mask heavily
    // while its peak still sits on the envelope, so the peak decides the status
    // and the overlap fraction is kept alongside as evidence.
    int[] shellCounts = shellCounts(nucleusOnly, peak, SHELL_RADIUS_UM, sx, sz, imp)
    double shellFraction = shellCounts[2] == 0 ? Double.NaN : (double) shellCounts[0] / shellCounts[2]
    // The bounds bracket every value the true fraction could take once the part
    // of the sphere that left the stack is admitted as unknown: all of it
    // cytoplasm, or all of it nucleus. shell_lower is the reported fraction,
    // because counting an unimaged direction as "not nucleus" is exactly the
    // lower bound; shell_upper says how far the truth could be above it. A row
    // whose bounds straddle ENGULFED_SHELL_FRACTION has an undetermined verdict,
    // not merely an uncertain one, and carries SHELL_TRUNCATED to say so.
    double shellOutOfVolume = 1.0 - (double) shellCounts[1] / shellCounts[2]
    double shellLower = (double) shellCounts[0] / shellCounts[2]
    double shellUpper = (double) (shellCounts[0] + shellCounts[2] - shellCounts[1]) / shellCounts[2]
    boolean shellTruncated = shellOutOfVolume > SHELL_OUT_OF_VOLUME_FRACTION
    double strictThreshold = autoThreshold(nucleiRaw, AutoThresholder.Method.Moments)
    double shellFractionStrict = shellFractionAboveThreshold(nucleiRaw, peak, SHELL_RADIUS_UM, sx, sz, strictThreshold)
    boolean thresholdSensitive = Math.abs(shellFraction - shellFractionStrict) > THRESHOLD_SENSITIVE_DELTA

    if (shellFraction >= ENGULFED_SHELL_FRACTION) {
        qcStatus = "ENGULFED_NOT_MEASURABLE"
    } else if (shellFraction >= DOCKED_SHELL_FRACTION) {
        qcStatus = "AT_NUCLEAR_ENVELOPE"
    } else {
        qcStatus = "OK"
    }
    if (thresholdSensitive) {
        qcStatus = qcStatus + "+THRESHOLD_SENSITIVE"
    }
    if (brighterIntranuclearRejected) {
        qcStatus = qcStatus + "+BRIGHTER_INTRANUCLEAR_REJECTED"
    }
    if (Math.abs(spotCenterToNucleusBorder) > IMPLAUSIBLE_DISTANCE_UM) {
        qcStatus = qcStatus + "+IMPLAUSIBLY_FAR"
    }
    // A centrosome cannot be inside the nucleus, so a negative distance is not a
    // measurement of anything: it says the spot centroid fell on the wrong side
    // of where the segmentation put the envelope. Every negative value in this
    // dataset is below 0.24 um, and both scales that govern the boundary are
    // larger: moving the nuclear threshold from Otsu to Moments shifts the
    // surface by 0.26 um on average (0.09-0.50 over the affected images), and the
    // axial resolution of the acquisition is about 0.5 um. So the sign carries no
    // information here; it is the segmentation's own uncertainty read out as a
    // number. The row is kept with its measured value rather than clamped, so
    // nothing is invented, and the flag lets the statistics drop it.
    if (spotCenterToNucleusBorder < 0) {
        qcStatus = qcStatus + "+NEGATIVE_DISTANCE"
    }
    // The measured surface was the crop, and the number reported is the one
    // recomputed against the real envelope.
    if (surfaceClipped) {
        qcStatus = qcStatus + "+SURFACE_CLIPPED"
    }
    // Enough of the envelope is missing that even an uncontaminated distance
    // describes a nucleus the stack only partly contains.
    if (nucleusTruncated) {
        qcStatus = qcStatus + "+NUCLEUS_TRUNCATED"
    }
    // shell_in_nucleus was renormalised over the part of the sphere that was
    // imaged; read shell_lower and shell_upper instead of the point estimate.
    if (shellTruncated) {
        qcStatus = qcStatus + "+SHELL_TRUNCATED"
    }
    // The runner-up matters: a second punctum a micron away is the other half of
    // a duplicated pair, and which half got measured is then arbitrary.
    def otherSpots = spotCandidates.findAll { !it.is(spot) }
    double secondSpotDistance = Double.NaN
    double pairSeparation = Double.NaN
    if (!otherSpots.isEmpty()) {
        Object3D runnerUp = otherSpots.min { it.distCenterUnit(spot) }
        pairSeparation = runnerUp.distCenterUnit(spot)
        secondSpotDistance = runnerUp.distCenterBorderUnit(nucleus)
        if (pairSeparation <= PAIR_SEPARATION_UM) {
            qcStatus = qcStatus + "+CENTROSOME_PAIR"
        }
    }
    if (spot.getVolumePixels() > PLAUSIBLE_SPOT_VOXELS) {
        qcStatus = qcStatus + "+OVERSIZED_SPOT"
    }
    if (nucleusCandidates.size() > 1) {
        qcStatus = qcStatus + "+MULTIPLE_NUCLEI"
    }
    // Written whenever the crop touched the measurement, so the run can be
    // reviewed by grepping the log instead of reopening every image.
    if (surfaceClipped || nucleusTruncated || shellTruncated) {
        logFile << "${file.name}: truncation ${String.format('%.1f', 100 * truncatedSurfaceFraction)}% of the nuclear surface is the edge of the stack; "
        logFile << "centroid ${String.format('%.4f', centroidToAll)} -> ${String.format('%.4f', centroidToFree)} um, "
        logFile << "peak ${String.format('%.4f', peakToAll)} -> ${String.format('%.4f', peakToFree)} um, "
        logFile << "border ${String.format('%.4f', borderToAll)} -> ${String.format('%.4f', borderToFree)} um; "
        logFile << "shell ${String.format('%.2f', shellFraction)} over ${String.format('%.0f', 100 * (1 - shellOutOfVolume))}% of the sphere, true value in [${String.format('%.2f', shellLower)}, ${String.format('%.2f', shellUpper)}]\n"
    }

    File csv = new File(groupDir, stem + "_spot_to_nucleus_distance.xls")
    csv.text = "file_name\tdistance_to_nucleus_surface_um\tnegative_distance\tdistance_from_peak_um\tsurface_to_surface_um\tqc_status\tshell_in_nucleus\tshell_strict\tdapi_at_spot_ratio\tspot_voxels\tspot_n\tnuclei_n\tpair_separation_um\tsecond_spot_um\tspot_in_nucleus_pct\tspot_z_clipped\tnucleus_z_clipped\tnucleus_surface_clipped_frac\tshell_out_of_volume\tshell_lower\tshell_upper\tspot_id\tnucleus_id\tspot_x_um\tspot_y_um\tspot_z_um\n"
    String row = [file.name, spotCenterToNucleusBorder, spotCenterToNucleusBorder < 0, peakToNucleusBorder, surfaceToSurface, qcStatus, String.format("%.2f", shellFraction), String.format("%.2f", shellFractionStrict), String.format("%.2f", dapiRatio), spot.getVolumePixels(), spotCandidates.size(), nucleusCandidates.size(),
            Double.isNaN(pairSeparation) ? "" : String.format("%.3f", pairSeparation),
            Double.isNaN(secondSpotDistance) ? "" : String.format("%.3f", secondSpotDistance), String.format("%.1f", overlapPercent), spotZClipped, nucleusZClipped,
            String.format("%.3f", truncatedSurfaceFraction), String.format("%.3f", shellOutOfVolume), String.format("%.2f", shellLower), String.format("%.2f", shellUpper),
            1, 1, spot.getCenterX() * sx, spot.getCenterY() * sy, spot.getCenterZ() * sz].join("\t") + "\n"
    csv << row
    summaryFile << row

    saveQcOverlay(imp, nucleusDisplay, spotDisplay, nucleiList, nucleus, spot, peak,
            "${stem}  ${qcStatus}  d=${String.format('%+.3f', spotCenterToNucleusBorder)}um",
            new File(overlayDir, "qc_" + stem + ".png"), sx, sz)

    // Label images are retained for visual inspection in Fiji.
    saveLabels(nucleiList, imp, new File(groupDir, stem + "_nuclei_labels.tif"), sx, sz)
    saveLabels(spotsList, imp, new File(groupDir, stem + "_spots_labels.tif"), sx, sz)
    new File(groupDir, stem + "_metadata.txt").text = "file=${file.name}\nchannels=${imp.getNChannels()}\nslices=${imp.getNSlices()}\nframes=${imp.getNFrames()}\nwidth=${imp.getWidth()}\nheight=${imp.getHeight()}\npixel_width_um=${sx}\npixel_height_um=${sy}\npixel_depth_um=${sz}\nspot_threshold=${spotThreshold}\nmin_spot_voxels=${MIN_SPOT_VOXELS}\nnucleus_candidates=${nucleusCandidates.size()}\nspot_candidates=${spotCandidates.size()}\nnucleus_threshold_method=${NUCLEUS_THRESHOLD_METHOD}\nspot_threshold_method=${SPOT_THRESHOLD_METHOD}\nmax_spot_voxels=${MAX_SPOT_VOXELS}\noversized_spots_rejected=${oversizedSpots.size()}\nspot_select_by=${SPOT_SELECT_BY}\nselected_spot_volume_voxels=${spot.getVolumePixels()}\nselected_spot_max_intensity=${spot.getPixMaxValue(spotIntensity)}\nqc_status=${qcStatus}\nnegative_distance=${spotCenterToNucleusBorder < 0}\naxial_resolution_um=${AXIAL_RESOLUTION_UM}\ndapi_at_spot=${dapiAtSpot}\ndapi_nucleus_median=${nucleusMedian}\ndapi_at_spot_ratio=${dapiRatio}\nspot_candidates_extranuclear=${extranuclear.size()}\nbrighter_intranuclear_rejected=${brighterIntranuclearRejected}\npair_separation_um=${pairSeparation}\nsecond_spot_distance_um=${secondSpotDistance}\nplausible_spot_voxels=${PLAUSIBLE_SPOT_VOXELS}\nshell_radius_um=${SHELL_RADIUS_UM}\nshell_in_nucleus=${shellFraction}\nshell_in_nucleus_strict=${shellFractionStrict}\nstrict_threshold=${strictThreshold}\nthreshold_sensitive=${thresholdSensitive}\ncentroid_inside_nucleus=${centroidInside}\npeak_inside_nucleus=${peakInside}\nspot_peak_voxel=${(int) peak.getX()},${(int) peak.getY()},${(int) peak.getZ()}\nspot_in_nucleus_pct=${overlapPercent}\nspot_overlap_voxels=${overlapVoxels}\nspot_z_range=${spot.getZmin()}..${spot.getZmax()}\nnucleus_z_range=${nucleus.getZmin()}..${nucleus.getZmax()}\nlast_slice_index=${lastSlice}\nnucleus_surface_voxels=${nucleusSurface.size()}\nnucleus_surface_clipped_voxels=${nucleusSurface.size() - nucleusSurfaceFree.size()}\nnucleus_surface_clipped_frac=${truncatedSurfaceFraction}\nnucleus_truncated=${nucleusTruncated}\ndist_centroid_to_any_surface_um=${centroidToAll}\ndist_centroid_to_free_surface_um=${centroidToFree}\ndist_peak_to_any_surface_um=${peakToAll}\ndist_peak_to_free_surface_um=${peakToFree}\ndist_border_to_any_surface_um=${borderToAll}\ndist_border_to_free_surface_um=${borderToFree}\nsurface_clipped=${surfaceClipped}\nshell_samples_total=${shellCounts[2]}\nshell_samples_in_volume=${shellCounts[1]}\nshell_out_of_volume=${shellOutOfVolume}\nshell_lower=${shellLower}\nshell_upper=${shellUpper}\nshell_truncated=${shellTruncated}\nselected_nucleus_count=${nucleiList.size()}\nselected_spot_count=${spotsList.size()}\n"

    nucleusOnly.closeImagePlus()
    spotOnly.closeImagePlus()
    nucleusDisplay.close()
    spotDisplay.close()
    nucleiRaw.close()
    spotsRaw.close()
    nucleiMask.close()
    spotsMask.close()
    imp.close()
    logFile << "${file.name}: completed\n"
    } catch (Throwable error) {
        String message = error.toString() + "\n" + error.stackTrace.join("\n")
        new File(new File(outputDir, groupOf(stem)), stem + "_ERROR.txt").text = message
        logFile << "${file.name}: ERROR ${error}\n"
        IJ.log("ERROR in ${file.name}: ${error}")
    }
}

IJ.log("Finished. Results: " + outputDir.absolutePath)

// Reproduces the 3D Manager "3D Distance" table for an ordered object list.
// One row per ordered pair, so two objects always give exactly two rows.
List<String> distanceRows(List<Object3D> objects, Objects3DPopulation population) {
    List<String> rows = []
    int nb = 0
    for (int i = 0; i < objects.size(); i++) {
        for (int j = 0; j < objects.size(); j++) {
            if (i == j) {
                continue
            }
            Object3D a = objects[i]
            Object3D b = objects[j]
            double cenCen = a.distCenterUnit(b)
            double cenBor = a.distCenterBorderUnit(b)
            double borBor = a.distBorderUnit(b)
            double radiusCen = a.radiusCenter(b)
            double radiusCenBack = b.radiusCenter(a)
            double borRad = cenCen - radiusCen - radiusCenBack
            double excen = radiusCen != 0 ? cenCen / radiusCen : Double.NaN
            double periph = radiusCen != 0 ? borRad / radiusCen : Double.NaN
            Object3D closestCen = population.closestCenter(a, true)
            Object3D closestBor = population.closestBorder(a)
            rows << ([nb, a.getValue(), b.getValue(), a.getType(), b.getType(), a.getName(), b.getName(),
                      cenCen, cenBor, borBor, radiusCen, excen, borRad, periph,
                      closestCen == null ? "" : population.getIndexOf(closestCen) + 1,
                      closestBor == null ? "" : population.getIndexOf(closestBor) + 1,
                      closestCen == null ? "" : closestCen.getName(),
                      closestBor == null ? "" : closestBor.getName()].join("\t") + "\n")
            nb++
        }
    }
    return rows
}

ImagePlus extractChannel(ImagePlus source, int channel) {
    ImageStack sourceStack = source.getStack()
    int width = source.getWidth()
    int height = source.getHeight()
    int slices = source.getNSlices()
    int frames = source.getNFrames()
    ImageStack result = new ImageStack(width, height)
    for (int t = 1; t <= frames; t++) {
        for (int z = 1; z <= slices; z++) {
            int index = source.getStackIndex(channel, z, t)
            result.addSlice(sourceStack.getProcessor(index).duplicate())
        }
    }
    ImagePlus output = new ImagePlus("C" + channel, result)
    output.setDimensions(1, slices, frames)
    output.setCalibration(source.getCalibration().copy())
    return output
}

ImagePlus thresholdStack(ImagePlus source, AutoThresholder.Method method, double fixedThreshold) {
    ImageStack result = new ImageStack(source.getWidth(), source.getHeight())
    StackStatistics statistics = new StackStatistics(source)
    int thresholdBin = method == null ? -1 : new AutoThresholder().getThreshold(method, intHistogram(statistics.histogram))
    double threshold = method == null ? fixedThreshold : statistics.histMin + (thresholdBin + 0.5) * statistics.binSize
    for (int index = 1; index <= source.getStackSize(); index++) {
        ImageProcessor ip = source.getStack().getProcessor(index)
        ByteProcessor mask = new ByteProcessor(ip.getWidth(), ip.getHeight())
        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                mask.set(x, y, ip.getf(x, y) >= threshold ? 255 : 0)
            }
        }
        result.addSlice(mask)
    }
    ImagePlus output = new ImagePlus(source.getTitle() + " mask", result)
    output.setDimensions(source.getNChannels(), source.getNSlices(), source.getNFrames())
    output.setCalibration(source.getCalibration().copy())
    return output
}

int[] intHistogram(long[] histogram) {
    int[] converted = new int[histogram.length]
    for (int index = 0; index < histogram.length; index++) {
        converted[index] = (int)Math.min(Integer.MAX_VALUE, histogram[index])
    }
    return converted
}

void saveLabels(List<Object3D> objects, ImagePlus reference, File target, double sx, double sz) {
    ImageStack labels = new ImageStack(reference.getWidth(), reference.getHeight())
    for (int index = 1; index <= reference.getNSlices() * reference.getNFrames(); index++) {
        labels.addSlice(new ShortProcessor(reference.getWidth(), reference.getHeight()))
    }
    ImagePlus labelImage = new ImagePlus(target.name, labels)
    labelImage.setDimensions(1, reference.getNSlices(), reference.getNFrames())
    labelImage.setCalibration(reference.getCalibration().copy())
    ImageHandler handler = ImageHandler.wrap(labelImage)
    objects.each { Object3D object -> object.draw(handler, object.getValue()) }
    new FileSaver(labelImage).saveAsTiff(target.absolutePath)
    labelImage.close()
}

// An empty 16-bit stack matching the reference geometry, used to rasterise one
// object on its own so that overlap can be counted voxel by voxel.
ImagePlus blankLike(ImagePlus reference) {
    ImageStack stack = new ImageStack(reference.getWidth(), reference.getHeight())
    for (int index = 0; index < reference.getNSlices(); index++) {
        stack.addSlice(new ShortProcessor(reference.getWidth(), reference.getHeight()))
    }
    ImagePlus output = new ImagePlus("blank", stack)
    output.setDimensions(1, reference.getNSlices(), 1)
    output.setCalibration(reference.getCalibration().copy())
    return output
}

// Tests one voxel of a rasterised object, treating anything outside the image
// as background rather than letting the coordinate clamp to an edge voxel.
boolean insideMask(ImageHandler mask, int x, int y, int z, ImagePlus reference) {
    if (x < 0 || x >= reference.getWidth() || y < 0 || y >= reference.getHeight() || z < 0 || z >= reference.getNSlices()) {
        return false
    }
    return mask.getPixel(x, y, z) > 0
}

// The surface of a rasterised object: every voxel with at least one 6-neighbour
// that is background, with everything past the edge of the image counted as
// background. That is the contour mcib3d's distCenterBorderUnit and
// distBorderUnit search, reproduced here so the same set can be filtered before
// the search. Checked against the saved label images of all 81 rows of the
// previous run: the nearest-contour distance matches mcib3d to 1e-10 um.
List<Voxel3D> surfaceVoxels(Object3D object, ImageHandler mask, ImagePlus reference) {
    List<Voxel3D> surface = []
    object.getVoxels().each { Voxel3D voxel ->
        int x = (int) voxel.getX(), y = (int) voxel.getY(), z = (int) voxel.getZ()
        if (!insideMask(mask, x, y, z, reference)) {
            return
        }
        if (!insideMask(mask, x - 1, y, z, reference) || !insideMask(mask, x + 1, y, z, reference)
                || !insideMask(mask, x, y - 1, z, reference) || !insideMask(mask, x, y + 1, z, reference)
                || !insideMask(mask, x, y, z - 1, reference) || !insideMask(mask, x, y, z + 1, reference)) {
            surface << voxel
        }
    }
    return surface
}

// Whether a surface voxel is surface only because the acquisition stopped: it
// lies on the first or last slice or on the XY border, and every neighbour that
// was actually imaged is foreground. The object continues past it, so the flat
// plane through voxels like this is the edge of the stack, not the envelope.
// A voxel on the border that also has a genuine background neighbour inside the
// image is a real piece of surface and is deliberately kept.
boolean isTruncationFace(Voxel3D voxel, ImageHandler mask, ImagePlus reference) {
    int x = (int) voxel.getX(), y = (int) voxel.getY(), z = (int) voxel.getZ()
    int w = reference.getWidth(), h = reference.getHeight(), d = reference.getNSlices()
    if (!(x == 0 || y == 0 || z == 0 || x == w - 1 || y == h - 1 || z == d - 1)) {
        return false
    }
    int[][] steps = [[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]]
    for (int[] step : steps) {
        int nx = x + step[0], ny = y + step[1], nz = z + step[2]
        boolean outsideImage = nx < 0 || ny < 0 || nz < 0 || nx >= w || ny >= h || nz >= d
        if (!outsideImage && mask.getPixel(nx, ny, nz) <= 0) {
            return false
        }
    }
    return true
}

// Calibrated distance from a point to the nearest voxel of a surface set.
double nearestSurfaceUnit(List<Voxel3D> surface, double x, double y, double z, double sx, double sy, double sz) {
    double best = Double.NaN
    surface.each { Voxel3D voxel ->
        double dx = (voxel.getX() - x) * sx
        double dy = (voxel.getY() - y) * sy
        double dz = (voxel.getZ() - z) * sz
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz)
        if (Double.isNaN(best) || distance < best) {
            best = distance
        }
    }
    return best
}

// Closest approach between two surface sets, the surface-to-surface distance.
double nearestSurfacePair(List<Voxel3D> from, List<Voxel3D> to, double sx, double sy, double sz) {
    double best = Double.NaN
    from.each { Voxel3D voxel ->
        double distance = nearestSurfaceUnit(to, voxel.getX(), voxel.getY(), voxel.getZ(), sx, sy, sz)
        if (!Double.isNaN(distance) && (Double.isNaN(best) || distance < best)) {
            best = distance
        }
    }
    return best
}

// Evenly spaced directions on a sphere via the Fibonacci lattice, so the shell
// is sampled without the pole clustering that latitude/longitude steps give.
List<double[]> sphereDirections(int count) {
    List<double[]> directions = []
    double golden = Math.PI * (1 + Math.sqrt(5.0))
    for (int i = 0; i < count; i++) {
        double phi = Math.acos(1 - 2 * ((i + 0.5) / count))
        double theta = golden * i
        directions << ([Math.sin(phi) * Math.cos(theta), Math.sin(phi) * Math.sin(theta), Math.cos(phi)] as double[])
    }
    return directions
}

// [samples in the mask, samples inside the acquired volume, samples cast]. The
// three counts are kept apart because a direction that leaves the stack is not
// a miss: nothing was imaged there, so it can be neither counted nor dismissed.
int[] shellCounts(ImageHandler mask, Voxel3D centre, double radius, double sx, double sz, ImagePlus reference) {
    int hits = 0, inVolume = 0, total = 0
    sphereDirections(400).each { double[] u ->
        total++
        int ix = (int) Math.round(centre.getX() + u[0] * radius / sx)
        int iy = (int) Math.round(centre.getY() + u[1] * radius / sx)
        int iz = (int) Math.round(centre.getZ() + u[2] * radius / sz)
        if (ix >= 0 && ix < reference.getWidth() && iy >= 0 && iy < reference.getHeight() && iz >= 0 && iz < reference.getNSlices()) {
            inVolume++
            if (mask.getPixel(ix, iy, iz) > 0) hits++
        }
    }
    return [hits, inVolume, total] as int[]
}

// A direction that leaves the stack counts against the fraction. Renormalising
// over only the imaged part of the sphere instead would divide by a shrinking
// denominator exactly where the missing part is the cytoplasmic side, which read
// a docked spot as engulfed: one crop came out at 0.90 over 64% of the sphere and
// at 0.58 over all of it. Nothing was recorded outside the stack, so this is not
// a measurement of what is there; it is the choice to let an unimaged direction
// count as "not nucleus", which is what makes the number comparable between a
// spot in the middle of the stack and one near its edge. shell_out_of_volume and
// shell_upper carry how much of the sphere was assumed rather than seen.
double shellFractionInMask(ImageHandler mask, Voxel3D centre, double radius, double sx, double sz, ImagePlus reference) {
    int[] counts = shellCounts(mask, centre, radius, sx, sz, reference)
    return counts[2] == 0 ? Double.NaN : (double) counts[0] / counts[2]
}

// The same convention as shellFractionInMask, so that the two fractions stay
// comparable and THRESHOLD_SENSITIVE keeps measuring threshold sensitivity
// rather than a difference in how the stack edge is treated.
double shellFractionAboveThreshold(ImagePlus source, Voxel3D centre, double radius, double sx, double sz, double threshold) {
    int hits = 0, total = 0
    sphereDirections(400).each { double[] u ->
        total++
        int ix = (int) Math.round(centre.getX() + u[0] * radius / sx)
        int iy = (int) Math.round(centre.getY() + u[1] * radius / sx)
        int iz = (int) Math.round(centre.getZ() + u[2] * radius / sz)
        if (ix >= 0 && ix < source.getWidth() && iy >= 0 && iy < source.getHeight() && iz >= 0 && iz < source.getNSlices()) {
            if (source.getStack().getProcessor(iz + 1).getf(ix, iy) >= threshold) hits++
        }
    }
    return total == 0 ? Double.NaN : (double) hits / total
}

double autoThreshold(ImagePlus source, AutoThresholder.Method method) {
    StackStatistics statistics = new StackStatistics(source)
    int bin = new AutoThresholder().getThreshold(method, intHistogram(statistics.histogram))
    return statistics.histMin + (bin + 0.5) * statistics.binSize
}

// Renders what was segmented: the spot peak seen in XY, XZ and YZ, with the
// nucleus kept outlined in magenta, any nucleus that was discarded in green,
// and the chosen spot in yellow, over blue DAPI and red gamma-tubulin. The two
// orthogonal sections matter more than the XY one here, because the axial
// direction is where a docked spot and an engulfed one are told apart.
void saveQcOverlay(ImagePlus reference, ImagePlus nucleusDisplay, ImagePlus spotDisplay,
                   List<Object3D> allNuclei, Object3D nucleus, Object3D spot, Voxel3D peak,
                   String caption, File target, double sx, double sz) {
    int W = reference.getWidth(), H = reference.getHeight(), D = reference.getNSlices()
    int px = (int) peak.getX(), py = (int) peak.getY(), pz = (int) peak.getZ()

    ImagePlus allImage = blankLike(reference)
    ImageHandler allMask = ImageHandler.wrap(allImage)
    allNuclei.each { it.draw(allMask, 255) }
    ImagePlus keptImage = blankLike(reference)
    ImageHandler keptMask = ImageHandler.wrap(keptImage)
    nucleus.draw(keptMask, 255)
    ImagePlus spotImage = blankLike(reference)
    ImageHandler spotMask = ImageHandler.wrap(spotImage)
    spot.draw(spotMask, 255)

    double[] blueRange = displayRange(nucleusDisplay)
    double[] redRange = displayRange(spotDisplay)

    ColorProcessor xy = overlayPanel(W, H, blueRange, redRange,
            { int a, int b -> [nucleusDisplay.getStack().getProcessor(pz + 1).getf(a, b), spotDisplay.getStack().getProcessor(pz + 1).getf(a, b)] as double[] },
            { int a, int b -> allMask.getPixel(a, b, pz) }, { int a, int b -> keptMask.getPixel(a, b, pz) }, { int a, int b -> spotMask.getPixel(a, b, pz) })
    // The XY view above sits at the spot's peak slice, which is where the spot
    // is but often near the top or bottom of the nucleus, where DAPI is dim and
    // the outline reads as too large. A second XY view through the middle of the
    // nucleus shows the mask where the nucleus is actually bright.
    int nz = Math.max(0, Math.min(D - 1, (int) Math.round(nucleus.getCenterZ())))
    ColorProcessor xyMid = overlayPanel(W, H, blueRange, redRange,
            { int a, int b -> [nucleusDisplay.getStack().getProcessor(nz + 1).getf(a, b), spotDisplay.getStack().getProcessor(nz + 1).getf(a, b)] as double[] },
            { int a, int b -> allMask.getPixel(a, b, nz) }, { int a, int b -> keptMask.getPixel(a, b, nz) }, { int a, int b -> spotMask.getPixel(a, b, nz) })
    ColorProcessor xz = overlayPanel(W, D, blueRange, redRange,
            { int a, int b -> [nucleusDisplay.getStack().getProcessor(b + 1).getf(a, py), spotDisplay.getStack().getProcessor(b + 1).getf(a, py)] as double[] },
            { int a, int b -> allMask.getPixel(a, py, b) }, { int a, int b -> keptMask.getPixel(a, py, b) }, { int a, int b -> spotMask.getPixel(a, py, b) })
    ColorProcessor yz = overlayPanel(H, D, blueRange, redRange,
            { int a, int b -> [nucleusDisplay.getStack().getProcessor(b + 1).getf(px, a), spotDisplay.getStack().getProcessor(b + 1).getf(px, a)] as double[] },
            { int a, int b -> allMask.getPixel(px, a, b) }, { int a, int b -> keptMask.getPixel(px, a, b) }, { int a, int b -> spotMask.getPixel(px, a, b) })
    drawCrosshair(xy, px, py); drawCrosshair(xz, px, pz); drawCrosshair(yz, py, pz)

    int zoom = 3
    // z is sampled finer than xy, so the sections are squeezed back to the same
    // micron-per-pixel as the XY view instead of appearing stretched.
    double zAspect = sz / sx
    def xyS = xy.resize(W * zoom, H * zoom, false)
    def xyMidS = xyMid.resize(W * zoom, H * zoom, false)
    def xzS = xz.resize(W * zoom, Math.max(1, (int) Math.round(D * zAspect * zoom)), false)
    def yzS = yz.resize(H * zoom, Math.max(1, (int) Math.round(D * zAspect * zoom)), false)

    int gap = 12, header = 20
    int leftHeight = xyS.getHeight() + gap + xyMidS.getHeight()
    int totalWidth = xyS.getWidth() + gap + Math.max(xzS.getWidth(), yzS.getWidth())
    int totalHeight = header + Math.max(leftHeight, xzS.getHeight() + gap + yzS.getHeight())
    ColorProcessor canvas = new ColorProcessor(totalWidth, totalHeight)
    canvas.setColor(java.awt.Color.BLACK)
    canvas.fill()
    canvas.insert(xyS, 0, header)
    canvas.insert(xyMidS, 0, header + xyS.getHeight() + gap)
    canvas.insert(xzS, xyS.getWidth() + gap, header)
    canvas.insert(yzS, xyS.getWidth() + gap, header + xzS.getHeight() + gap)
    canvas.setColor(java.awt.Color.WHITE)
    canvas.drawString(caption, 4, 15)
    canvas.drawString("XY z=" + pz + " (spot peak)", 4, header + 15)
    canvas.drawString("XY z=" + nz + " (nucleus mid)", 4, header + xyS.getHeight() + gap + 15)
    canvas.drawString("XZ y=" + py, xyS.getWidth() + gap + 4, header + 15)
    canvas.drawString("YZ x=" + px, xyS.getWidth() + gap + 4, header + xzS.getHeight() + gap + 15)
    new FileSaver(new ImagePlus("qc", canvas)).saveAsPng(target.absolutePath)

    allMask.closeImagePlus(); keptMask.closeImagePlus(); spotMask.closeImagePlus()
}

// Script-level constants are locals in Groovy and are not visible inside a
// method, so the display gamma lives here rather than beside the segmentation
// parameters at the top. It affects the QC rendering only.
double gammaScale(double value, double[] range) {
    final double DISPLAY_GAMMA = 0.5
    double t = (value - range[0]) / (range[1] - range[0])
    if (t <= 0) return 0.0
    if (t >= 1) return 1.0
    return Math.pow(t, DISPLAY_GAMMA)
}

ColorProcessor overlayPanel(int w, int h, double[] blueRange, double[] redRange,
                            Closure intensities, Closure allMask, Closure keptMask, Closure spotMask) {
    ColorProcessor panel = new ColorProcessor(w, h)
    for (int b = 0; b < h; b++) {
        for (int a = 0; a < w; a++) {
            double[] v = intensities(a, b)
            // Gamma, not a wider clip. Segmentation runs on raw values with one
            // global threshold, so a slice near the top of the nucleus is mostly
            // dim voxels that still sit above the cut; displayed linearly they
            // read as black and the outline looks drawn around nothing (seen in one
            // at z=43: half the masked pixels render under 10% brightness).
            // Gamma lifts those without blowing out the bright core.
            int blue = toByte(gammaScale(v[0], blueRange) * 255)
            int red = toByte(gammaScale(v[1], redRange) * 255)
            panel.set(a, b, (red << 16) | blue)
        }
    }
    for (int b = 0; b < h; b++) {
        for (int a = 0; a < w; a++) {
            if (isContour(allMask, a, b, w, h)) panel.set(a, b, 0x00C000)
            if (isContour(keptMask, a, b, w, h)) panel.set(a, b, 0xFF00FF)
            if (isContour(spotMask, a, b, w, h)) panel.set(a, b, 0xFFFF00)
        }
    }
    return panel
}

boolean isContour(Closure mask, int a, int b, int w, int h) {
    if (mask(a, b) <= 0) return false
    if (a == 0 || b == 0 || a == w - 1 || b == h - 1) return true
    return mask(a - 1, b) <= 0 || mask(a + 1, b) <= 0 || mask(a, b - 1) <= 0 || mask(a, b + 1) <= 0
}

void drawCrosshair(ColorProcessor panel, int a, int b) {
    panel.setColor(java.awt.Color.WHITE)
    panel.drawLine(a - 3, b, a - 1, b); panel.drawLine(a + 1, b, a + 3, b)
    panel.drawLine(a, b - 3, a, b - 1); panel.drawLine(a, b + 1, a, b + 3)
}

int toByte(double value) {
    return (int) Math.max(0, Math.min(255, Math.round(value)))
}

// Clipping the top 0.05% keeps one saturated punctum from crushing everything
// else to black.
double[] displayRange(ImagePlus source) {
    StackStatistics statistics = new StackStatistics(source)
    long total = 0
    for (long count : statistics.histogram) total += count
    long target = (long) (total * 0.9995)
    long running = 0
    int bin = statistics.histogram.length - 1
    for (int i = 0; i < statistics.histogram.length; i++) {
        running += statistics.histogram[i]
        if (running >= target) { bin = i; break }
    }
    double high = statistics.histMin + (bin + 1) * statistics.binSize
    return [statistics.min, Math.max(high, statistics.min + 1)] as double[]
}

// Sample group: the leading letters of the name, so every SampleA1-* and SampleA2-*
// file is filed under SampleA. A name that starts with no letters forms its own
// group rather than landing loose.
String groupOf(String stem) {
    def match = stem =~ /^([A-Za-z]+)/
    return match ? match[0][1] : stem
}

// Sorts SampleA5-9 before SampleA5-10 by zero-padding every run of digits, which plain
// string order would get wrong.
String naturalKey(String name) {
    return name.toLowerCase().replaceAll(/\d+/) { String digits -> digits.padLeft(12, '0') }
}

// The detector saturates at its full-scale value, so a punctum's brightest
// voxels form a plateau rather than a single maximum. Object3D.getPixelMax
// returns the first voxel reaching that maximum, and because the voxel list is
// built slice by slice that is always the plateau's lowest-z voxel: measured
// over this dataset it sat 4.15 slices (0.50 um) below the object centroid in
// 39 of 40 images, a bias large enough to move a spot across the envelope. The
// centre of the plateau is the stable anchor.
Voxel3D brightestPlateauCentre(Object3D object, ImageHandler intensity) {
    double maximum = object.getPixMaxValue(intensity)
    double cut = maximum * 0.99
    double sx = 0, sy = 0, sz = 0
    int count = 0
    object.getVoxels().each { Voxel3D voxel ->
        if (intensity.getPixel((int) voxel.getX(), (int) voxel.getY(), (int) voxel.getZ()) >= cut) {
            sx += voxel.getX(); sy += voxel.getY(); sz += voxel.getZ(); count++
        }
    }
    if (count == 0) {
        return object.getPixelMax(intensity)
    }
    return new Voxel3D(sx / count, sy / count, sz / count, maximum)
}
