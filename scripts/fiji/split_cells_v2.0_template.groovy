// PUBLIC TEMPLATE | crop workflow 2.0 | 2026-09-22
// Multichannel z-stack -> single-cell crops with selective foreign C1/C3 cleanup.
// This public template contains no experimental data and must be configured locally.
// 除版本文字、报告目录和 QC 标题外，算法、参数、人工保留规则及输出格式与 1.6 完全一致。
// 1.5.2 支持按人工目检把个别对象重新纳入（见下方 MANUAL_KEEP_*）；输出 tif 只裁 XY，
// 保留所有通道、所有 Z 层、16-bit 灰度和原标定，不做结构掩膜或像素置零。
// 1.5.3 复制 LUT，但源 ND2 的 Color 模式仍让 Fiji 单通道显示；1.5.4 复制 Info。
// 1.5.5 在 copyLuts 后再次强制 Composite，防止 copyLuts 把 Color 模式写回 TIFF。
// 1.6 对裁切框中明确属于其他核/其他中心体标签的 C1/C3 体素做最小清理；背景、目标结构、
// 其他通道和全 Z 均保留。显示 LUT/范围直接采用参考 TIFF，使观感一致。
//
// 1.4.1 的两项修复：
// (1) 递归分水岭只接受至少两个块都 >= 拆分前中位体积 30% 的切分。这样拒绝
//     那种 99%/1% 的伪切分，也阻止后续把大块切成 96%/4% 和 98%/1% 的碎片。
// (2) 1.5 不因邻核/邻中心体进入矩形裁切而剔除；输出 tif 时做结构掩膜，只保留目标核
//     和其归属的 γ-tubulin，其他结构置零，保证下游看到正好一核一中心体。邻近结构仍写 QC 标记。
//     仍疑似粘连的大核、无自身中心体、亮度不足或距离过远的对象继续 dropped。
//
// ============================ 与 1.3 的差别 ============================
//
// (1) 分水岭拆分改为递归（用户 2026-09-20 看过 5_split_qc.png 后决定，不降 tol）。
//
//     1.3 只拆一轮：5.nd2 的 v1.0 #8 是 3.79 x 中位体积，拆成 287787 / 54420 之后，
//     大块 287787 仍然是 2.87 x 中位体积，明显还是两个核黏在一起，但 1.3 到此为止。
//
//     1.4 改成：一轮拆完后逐块检查，只要某块体积仍然 > SPLIT_VOLUME_RATIO x 中位体积，
//     就对该块单独再跑一次 EDT + MaximaFinder + Watershed3D。
//       - 中位体积始终用“拆分前、已过 MIN_NUCLEUS_VOXELS 的对象”的那一个，不每轮重算，
//         这样判据在整个递归过程里是同一条线，不会因为拆出小块把中位数拉低而越拆越碎。
//       - 最大递归深度 SPLIT_MAX_DEPTH = 3。
//       - 终止条件：某块的 MaximaFinder 只找到 1 个种子，或分水岭没有把它分成 2 块以上，
//         就停止对该块继续拆（该块原样保留，包括它自己的 dam 体素）。
//       - 每一层的每一次实际拆分都记进 <stem>_split_objects.csv，带 depth 列。
//
//     块的命名：母对象的 v1.0 编号 + 每层一个后缀字符，按该层体积从大到小分配。
//       第 1 层用 a b c...，第 2 层用 1 2 3...，第 3 层用 x y z...
//       例如 8 -> 8a / 8b，8a 再拆 -> 8a1 / 8a2。
//
//     tol、半径、只对超标对象动手、拆完重新做体积筛选，全部与 1.3 相同，未改。
//
// (2) 目检图的排版修掉（用户 2026-09-20 指出 dropped/4_dropped.png 文字互相压住、读不了）。
//
//     a. 所有 drawString 一律先用 ColorProcessor.getStringWidth()（即 FontMetrics 实测）
//        量宽度，超过可用宽度就截断加 "..."，不再凭估计排版。截断次数写进日志，
//        以便发现“图上信息被悄悄截掉”。
//     b. 被剔除细胞缩略图的标注改成每项一行（编号 / 核体积 / 最亮候选 / 剔除理由），
//        不再把两项塞进同一行。格子最小宽度提到 DROPPED_TILE_MIN_W = 230。
//     c. 候选标记大幅减量：一个格子里只画一个空心圆环 —— 该核自己分到的最亮候选
//        （若一个都没分到，就画窗口内最亮的那个，并标出它归属给了谁）。其余候选只画
//        2 px 的小点，不带任何文字。1.3 在 #19 那格画了 22 个带标签的圆环，把细胞盖死了。
//     d. 归属编号只在那一个圆环旁标一次。
//     e. 表头说明改成两行。
//     f. kept/<stem>_overview.png 的每一行文字同样走实测宽度截断，格子最小宽度提到 200。
//
//     保持不变：青色描出被剔除核自身的 XY 投影轮廓；格子上仍标旧编号、核体积、EDGE；
//     圆环仍然画在信号外侧，不盖住红色中心体信号。
//
// 其余一律沿用 1.3：
//   - 保留判据 = 该核分到的候选中最亮的那个满足 maxI >= MIN_SPOT_MAX_INTENSITY(2000)
//     且到本核表面距离 <= MAX_OWN_SPOT_DISTANCE_UM(5.0 um)。2000 取自 4.nd2 实测空档
//     （保留组最亮候选最低 3559，剔除组最亮候选最高 1434）的中部。
//   - 核分割参数 sigma=1.0/1.0、全栈 Otsu、MIN_NUCLEUS_VOXELS=1000、CROP_MARGIN_PX=12、
//     ROW_BAND_PX=60；逐核局部 MaxEntropy 检出；并集掩膜后做一次全图连通域去重。
//   - 输出目录 <OUTPUT_DIR>/kept、<OUTPUT_DIR>/dropped、<OUTPUT_DIR>/qc。
//   - CSV 与日志写 archive/20260920_split_v2.0/，不覆盖旧版表。
//
// 本脚本只做分割与裁切，不产生任何距离测量结果：距离仍然只由
// batch_3d_suite_distance.groovy 负责，且要等用户目检完再跑，本脚本不代劳。
// 不读写结果目录，也绝不写回正式输入目录或合并数据目录。
//
// 运行：
//   fiji --headless --console --run split_cells_v2.0_template.groovy
import ij.IJ
import ij.ImagePlus
import ij.CompositeImage
import ij.ImageStack
import ij.measure.Calibration
import ij.plugin.GaussianBlur3D
import ij.process.ByteProcessor
import ij.process.ColorProcessor
import ij.process.ImageProcessor
import ij.process.StackStatistics
import ij.process.AutoThresholder
import ij.io.FileSaver
import loci.plugins.BF
import mcib3d.image3d.ImageHandler
import mcib3d.image3d.ImageInt
import mcib3d.image3d.ImageByte
import mcib3d.image3d.ImageShort
import mcib3d.image3d.ImageLabeller
import mcib3d.image3d.distanceMap3d.EDT
import mcib3d.image3d.processing.MaximaFinder
import mcib3d.image3d.regionGrowing.Watershed3D
import mcib3d.geom.Object3D
import mcib3d.geom.Voxel3D
import mcib3d.geom.Objects3DPopulation
import java.awt.Color
import java.awt.Font

// ===================== 可调参数与输入输出路径 =====================

def projectDir = new File(System.getProperty("projectDir", "."))
// 本次要处理的原图。换细胞系时改这一行和下面三行路径即可。
def INPUT_FILES = [new File(projectDir, "data/combined_data/sample_01.nd2")]
def DISPLAY_REFERENCE_FILE = new File(projectDir, "data/reference/display_reference.tif")
def SPLIT_DIR = new File(projectDir, "data/split_data/split_cells")
def OUTPUT_DIR = new File(SPLIT_DIR, "kept")        // 单细胞 tif + 总览图
def DROPPED_DIR = new File(SPLIT_DIR, "dropped")    // 被剔除细胞缩略图
def QC_DIR = new File(SPLIT_DIR, "qc")              // 分水岭拆分前后对照图
def REPORT_DIR = new File(projectDir, "results/split_report_v2.0") // CSV + logs
final String NAME_PREFIX = "CELL"

// 核通道（C1 = DAPI）与中心体通道（C3 = γ-tubulin），与测距流程相同。
final int NUCLEUS_CHANNEL = 1
final int SPOT_CHANNEL = 3

// ---- 核分割：以下八个常量与 1.0 / 1.1 / 1.2 逐字相同，本版未改动 ----
final double NUCLEUS_SIGMA_XY = 1.0
final double NUCLEUS_SIGMA_Z = 1.0
final AutoThresholder.Method NUCLEUS_THRESHOLD_METHOD = AutoThresholder.Method.Otsu
final int MIN_NUCLEUS_VOXELS = 1000
final int REPORT_MIN_VOXELS = 100
final double FUSED_VOLUME_RATIO = 1.8
final int CROP_MARGIN_PX = 12
final double ROW_BAND_PX = 60.0

// ---- 1.3 新增：粘连核 3D 分水岭拆分 ----
final boolean DO_WATERSHED_SPLIT = true
// 只有体积 > 该倍数 x 拆分前中位体积的对象才送去拆分，与 FUSED? 判据同一个倍数。
final double SPLIT_VOLUME_RATIO = FUSED_VOLUME_RATIO
// MaximaFinder 的种子搜索半径，像素单位。本批标定 XY=0.2859 um、Z=0.12 um，
// 2 px XY = 0.57 um、5 px Z = 0.60 um，两个方向的物理半径基本相同。
final float SPLIT_MAXIMA_RAD_XY = 2.0f
final float SPLIT_MAXIMA_RAD_Z = 5.0f
// 颈部判据：EDT 值 = 内切球半径，tol 要求颈部比两核半径都细该数值（um）以上才拆。
// 探针实测 0.5 / 1.0 过拆（3–11 块），2.0 一个都不拆，1.5 与目检一致。
final float SPLIT_NOISE_TOLERANCE_UM = 1.5f
// 分水岭只淹没 EDT > 该值的体素，即只在核内部蔓延，不进背景。
final double SPLIT_WATERSHED_THRESHOLD = 0.0001d
// 取子体积时四周留的空白，保证 EDT 在对象外侧确实是 0。
final int SPLIT_PAD_PX = 2
// 1.4 新增：递归拆分的最大层数。第 1 层是对原始连通域下刀，第 2 层是对第 1 层拆出的
// 超标块再下刀，依此类推。3 层足以把 3.79 x 中位体积的对象拆到每块都低于门槛。
final int SPLIT_MAX_DEPTH = 3
// 分水岭切分只有至少两块都达到拆分前中位体积的 30% 时才接受；否则保留母块并标为未解决粘连。
final double SPLIT_ACCEPT_MIN_FRACTION_OF_PREMEDIAN = 0.30

// ---- 中心体候选：以下五个常量取自冻结基线 1.0，不要自行改动 ----
// 1.2 把作用范围改回基线的语义：每个核的裁切窗口各算一次 MaxEntropy。
final double SPOT_SIGMA_XY = 0.5
final double SPOT_SIGMA_Z = 0.5
final AutoThresholder.Method SPOT_THRESHOLD_METHOD = AutoThresholder.Method.MaxEntropy
final int MIN_SPOT_VOXELS = 20
final int MAX_SPOT_VOXELS = 10000
// 保留条件里的距离上限。与基线 IMPLAUSIBLE_DISTANCE_UM 一致。
final double MAX_OWN_SPOT_DISTANCE_UM = 5.0
// ---- 1.3 新增：保留条件里的亮度下限 ----
// 4.nd2 实测：保留组最亮候选最低 3559，剔除组最亮候选最高 1434，中间 2125 灰阶空档里
// 没有任何候选。2000 取该空档中部，两侧各留 500+ 灰阶余量。12 bit 数据饱和值 4095。
final double MIN_SPOT_MAX_INTENSITY = 2000.0
// 单细胞裁切中保留邻核或邻中心体会使后续测距的核/点选择失去一对一含义，转入 dropped。
// 邻近结构只写 QC 标记，不作为单独 drop 条件；tif 输出为未修改的原始矩形裁切。
final boolean REQUIRE_CLEAN_SINGLE_CELL_CROP = false
final boolean MASK_NON_TARGET_STRUCTURES = false
// 用户于 2026-09-20 目检确认的人工保留：按原始 v1.0 母核编号记录。
final Set MANUAL_KEEP_SERIES_4 = [] as Set
final Set MANUAL_KEEP_SERIES_5 = [] as Set
final int MIN_FOREIGN_NUCLEUS_VOXELS_IN_CROP = 1000

// 总览图（保留细胞）
final boolean MAKE_OVERVIEW = true
final int OVERVIEW_COLUMNS = 5
final int OVERVIEW_PAD = 8
final int OVERVIEW_HEADER = 56
final int OVERVIEW_MIN_TILE_W = 200
final double OVERVIEW_CLIP_LOW = 0.005
final double OVERVIEW_CLIP_HIGH = 0.995

// 被剔除细胞缩略图
final boolean MAKE_DROPPED_PREVIEW = true
// 固定窗口边长：被剔除的核没有输出 tif，用固定窗口而不是它的外接框，这样各格
// 尺度一致、好横向比较。140 px = 40.0 µm，能装下一个典型细胞核加一圈邻域。
final int DROPPED_WINDOW_PX = 140
final int DROPPED_COLUMNS = 5
final int DROPPED_MARKER_RADIUS_PX = 9   // 候选圆环半径；中心体本身半径约 2–4 px
// 1.4：格子最小宽度，保证四行标注排得下；不够宽时文字会被实测截断而不是串格。
final int DROPPED_TILE_MIN_W = 230
// 没被画圆环的其余候选只画这么大的小点，不带文字。
final int DROPPED_DOT_RADIUS_PX = 2

// 拆分对照图
final boolean MAKE_SPLIT_QC = true
final int SPLIT_QC_MARGIN_PX = 20        // 母对象外接框外扩多少像素做窗口
// 原始整视野决策图：绿色细轮廓 = kept，黄色细轮廓 = dropped；不在细胞内部写字。
final boolean MAKE_DECISION_MAP = true

// 关掉它就只出 CSV 和图，不写 tif。
final boolean SAVE_TIFS = true

// ===================== 以下为实现 =====================

// 被实测宽度截断过的字符串条数；每张图画完都记进日志，免得信息被悄悄截掉而没人知道。
@groovy.transform.Field int TRUNCATED_STRINGS = 0

OUTPUT_DIR.mkdirs()
DROPPED_DIR.mkdirs()
QC_DIR.mkdirs()
REPORT_DIR.mkdirs()

ImagePlus displayReference = IJ.openImage(DISPLAY_REFERENCE_FILE.absolutePath)
if (displayReference == null || displayReference.getNChannels() != 4) {
    throw new IllegalStateException("Display reference is missing or not 4-channel: " + DISPLAY_REFERENCE_FILE.absolutePath)
}

INPUT_FILES.each { File INPUT_FILE ->

String stem = INPUT_FILE.name.replaceFirst(/(?i)\.(nd2|tif|tiff)$/, "")
String namePrefix = NAME_PREFIX + stem + "-"
def logFile = new File(REPORT_DIR, "split_log_" + stem + ".txt")
logFile.text = ""
def log = { String message ->
    IJ.log(message)
    logFile << message + "\n"
}

log("VERSION 2.0 (frozen v1.6 algorithm: manual keep + selective foreign C1/C3 cleanup + reference Composite display) | input : " + INPUT_FILE.absolutePath)
log("kept tif + overview -> " + OUTPUT_DIR.absolutePath)
log("dropped preview     -> " + DROPPED_DIR.absolutePath)
log("split qc            -> " + QC_DIR.absolutePath)
log("csv + log           -> " + REPORT_DIR.absolutePath)

// 1.4 的保留集合和编号会再变一次，旧 tif 留在目录里会和新 tif 混在一起。
// 只删本前缀的 tif，不碰目录里其他东西。备份见
// archive/20260920_split_v2.0/ 保留本版全部记录；本目录独立，不混用旧产物。
def stale = OUTPUT_DIR.listFiles().findAll {
    it.isFile() && it.name.startsWith(namePrefix) && it.name.toLowerCase().endsWith(".tif")
}
stale.each { it.delete() }
log("removed ${stale.size()} stale ${namePrefix}*.tif from a previous version")
if (!INPUT_FILE.exists()) {
    throw new IllegalArgumentException("Input not found: " + INPUT_FILE.absolutePath)
}

ImagePlus imp = INPUT_FILE.name.toLowerCase().endsWith(".nd2") ?
        BF.openImagePlus(INPUT_FILE.absolutePath)[0] : IJ.openImage(INPUT_FILE.absolutePath)
if (imp == null) {
    throw new IllegalStateException("Could not open " + INPUT_FILE.absolutePath)
}

int width = imp.getWidth()
int height = imp.getHeight()
int channels = imp.getNChannels()
int slices = imp.getNSlices()
int frames = imp.getNFrames()
int bitDepth = imp.getBitDepth()
Calibration cal = imp.getCalibration().copy()
double sx = cal.pixelWidth
double sy = cal.pixelHeight
double sz = cal.pixelDepth
double voxelVolume = sx * sy * sz
log("source: ${width}x${height} px, C=${channels}, Z=${slices}, T=${frames}, ${bitDepth}-bit")
log("calibration: ${sx} x ${sy} x ${sz} ${cal.getUnit()}  (voxel = ${voxelVolume} ${cal.getUnit()}^3)")
if (channels < SPOT_CHANNEL) {
    throw new IllegalArgumentException("Fewer than ${SPOT_CHANNEL} channels in ${INPUT_FILE.name}")
}
if (!(sx > 0) || !(sy > 0) || !(sz > 0)) {
    throw new IllegalArgumentException("Missing voxel size in " + INPUT_FILE.name)
}

// --- 1) 核通道 -> 3D 高斯模糊 -> 全栈 Otsu -> 3D 连通域（与 1.0/1.1/1.2 相同）---
ImagePlus nucleiRaw = extractChannel(imp, NUCLEUS_CHANNEL)
GaussianBlur3D.blur(nucleiRaw, NUCLEUS_SIGMA_XY, NUCLEUS_SIGMA_XY, NUCLEUS_SIGMA_Z)
double nucleusThreshold = autoThreshold(nucleiRaw, NUCLEUS_THRESHOLD_METHOD)
log("nucleus threshold (${NUCLEUS_THRESHOLD_METHOD}) = ${nucleusThreshold}")
ImagePlus nucleiMask = thresholdStack(nucleiRaw, nucleusThreshold)

ImageLabeller labeller = new ImageLabeller()
labeller.setMinSize(1)
ImageInt preLabels = labeller.getLabels(ImageHandler.wrap(nucleiMask), false)
Objects3DPopulation prePopulation = new Objects3DPopulation(preLabels, 0)
prePopulation.setCalibration(sx, sz, cal.getUnit())
def preObjects = prePopulation.getObjectsList()
log("connected components (>=1 voxel): ${preObjects.size()}")

// --- 2) 拆分前的体积分布记录（与 1.0/1.1/1.2 相同的 components_all 表）---
def reported = preObjects.findAll { it.getVolumePixels() >= REPORT_MIN_VOXELS }
                         .sort { -it.getVolumePixels() }
def componentsCsv = new File(REPORT_DIR, stem + "_components_all.csv")
componentsCsv.text = "rank,volume_voxels,volume_um3,bbox_x0,bbox_x1,bbox_y0,bbox_y1,bbox_z0,bbox_z1," +
        "bbox_w_px,bbox_h_px,centroid_x_px,centroid_y_px,centroid_z_px,kept\n"
reported.eachWithIndex { Object3D object, int index ->
    int volume = (int) object.getVolumePixels()
    componentsCsv << [index + 1, volume, String.format("%.3f", volume * voxelVolume),
                      object.getXmin(), object.getXmax(), object.getYmin(), object.getYmax(),
                      object.getZmin(), object.getZmax(),
                      object.getXmax() - object.getXmin() + 1, object.getYmax() - object.getYmin() + 1,
                      String.format("%.2f", object.getCenterX()), String.format("%.2f", object.getCenterY()),
                      String.format("%.2f", object.getCenterZ()),
                      volume >= MIN_NUCLEUS_VOXELS ? 1 : 0].join(",") + "\n"
}
log("components >= ${REPORT_MIN_VOXELS} voxels written to ${componentsCsv.name}: ${reported.size()}")

// 拆分前的体积筛选：只用来 (a) 定 v1.0 编号，(b) 定拆分的中位体积基准。
def preKept = preObjects.findAll { it.getVolumePixels() >= MIN_NUCLEUS_VOXELS }
if (preKept.isEmpty()) {
    imp.close()
    throw new IllegalStateException("No nucleus passed MIN_NUCLEUS_VOXELS=${MIN_NUCLEUS_VOXELS}")
}
double preMedianVolume = medianOf(preKept.collect { (double) it.getVolumePixels() })
log("pre-split nuclei passing volume (>= ${MIN_NUCLEUS_VOXELS}): ${preKept.size()}, " +
    "median volume ${(int) preMedianVolume} voxels")

// v1.0 编号在拆分前的集合上重建，拆分出来的块继承母对象编号并加后缀。
def orderedPre = rowMajorOrder(preKept, ROW_BAND_PX)
def parentOldIdOfPreLabel = [:]
orderedPre.eachWithIndex { Object3D object, int index -> parentOldIdOfPreLabel[object.getValue()] = index + 1 }
log("v1.0 numbering reconstructed for ${orderedPre.size()} pre-split nuclei")

// --- 3) 粘连核 3D 分水岭拆分（1.3 的改动 2）---
double splitVolumeCutoff = SPLIT_VOLUME_RATIO * preMedianVolume
def splitCandidates = preKept.findAll { it.getVolumePixels() > splitVolumeCutoff }
                             .sort { parentOldIdOfPreLabel[it.getValue()] }
log("watershed split: cutoff = ${SPLIT_VOLUME_RATIO} x median = ${(int) splitVolumeCutoff} voxels; " +
    "${splitCandidates.size()} object(s) above it -> v1.0 ids " +
    splitCandidates.collect { parentOldIdOfPreLabel[it.getValue()] })

ImageShort nucleusLabels = new ImageShort("nuclei labels", width, height, slices)
int nextLabel = 0
def parentOldIdOf = [:]     // new label -> v1.0 母对象编号（整数）
def oldIdTextOf = [:]       // new label -> "9" 或 "8a" 或 "8a1"
def splitRecords = []       // 每个进入拆分判定的母对象一条（含它整棵递归树）
def splitEvents = []        // 每一次真正发生的拆分一条，带 depth

preObjects.each { Object3D object ->
    int preLabel = object.getValue()
    boolean isCandidate = DO_WATERSHED_SPLIT && splitCandidates.any { it.getValue() == preLabel }
    if (!isCandidate) {
        nextLabel++
        object.draw(nucleusLabels, nextLabel)
        Integer pid = parentOldIdOfPreLabel[preLabel]
        parentOldIdOf[nextLabel] = pid
        oldIdTextOf[nextLabel] = pid == null ? "" : ("" + pid)
        return
    }

    int pid = parentOldIdOfPreLabel[preLabel]
    int rootVox = (int) object.getVolumePixels()

    // ---- 递归拆分：按层推进，levels[d] 是第 d 轮之后现存的所有块 ----
    // 第 0 层只有母对象本身；每一轮对该层中体积仍然超标的块各跑一次分水岭。
    // 中位体积用的始终是 preMedianVolume（拆分前的那个），判据在整个递归中不变。
    def levels = []
    levels << [[id: "" + pid, voxels: object.getVoxels(), depth: 0]]
    int totalDam = 0
    for (int d = 0; d < SPLIT_MAX_DEPTH; d++) {
        def current = levels[d]
        def next = []
        boolean anySplit = false
        current.each { piece ->
            if (piece.voxels.size() <= splitVolumeCutoff) { next << piece; return }
            def res = watershedVoxels(piece.voxels, SPLIT_PAD_PX, width, height, slices,
                                      (float) sx, (float) sz, SPLIT_MAXIMA_RAD_XY,
                                      SPLIT_MAXIMA_RAD_Z, SPLIT_NOISE_TOLERANCE_UM,
                                      SPLIT_WATERSHED_THRESHOLD)
            int acceptedLargeParts = res.parts.count { it.size() >= SPLIT_ACCEPT_MIN_FRACTION_OF_PREMEDIAN * preMedianVolume }
            if (res.parts.size() <= 1 || acceptedLargeParts < 2) {
                // 没拆开，或只切出边缘碎片：母块原样保留，不再递归吞掉正常核。
                piece.stopped = true
                piece.seeds = res.nSeeds
                next << piece
                log("  v1.0 #${pid} depth ${d + 1}: piece ${piece.id} " +
                    "${piece.voxels.size()} vox (${String.format('%.2f', piece.voxels.size() / preMedianVolume)} x median), " +
                    "${res.nSeeds} seed(s), ${res.parts.size()} part(s), ${acceptedLargeParts} acceptable " +
                    "(>= ${String.format('%.0f', SPLIT_ACCEPT_MIN_FRACTION_OF_PREMEDIAN * 100)}% pre-median) -> NOT SPLIT, stop here")
                return
            }
            anySplit = true
            totalDam += res.dam
            def children = []
            res.parts.eachWithIndex { vl, int k ->
                String cid = piece.id + suffixChar(d + 1, k)
                children << [id: cid, voxels: vl, depth: d + 1]
            }
            next.addAll(children)
            splitEvents << [oldId: pid, depth: d + 1, parentId: piece.id,
                            parentVox: piece.voxels.size(), rootVox: rootVox,
                            nSeeds: res.nSeeds, nParts: res.parts.size(), dam: res.dam,
                            parts: children.collect { [id: it.id, voxels: it.voxels.size()] }]
            log("  v1.0 #${pid} depth ${d + 1}: ${piece.id} ${piece.voxels.size()} vox " +
                "(${String.format('%.2f', piece.voxels.size() / preMedianVolume)} x median), " +
                "${res.nSeeds} seed(s) -> ${res.parts.size()} part(s) " +
                children.collect { it.id + '=' + it.voxels.size() } + " + ${res.dam} dam voxels")
        }
        levels << next
        if (!anySplit) break
    }

    def leaves = levels.last()
    // 叶子块按体积从大到小落标签，方便看表。
    leaves.sort { -it.voxels.size() }
    def partInfos = []
    leaves.each { piece ->
        nextLabel++
        piece.label = nextLabel
        parentOldIdOf[nextLabel] = pid
        oldIdTextOf[nextLabel] = piece.id
        partInfos << [label: nextLabel, oldIdText: piece.id, voxels: piece.voxels.size(),
                      depth: piece.depth]
        piece.voxels.each { Voxel3D v ->
            nucleusLabels.setPixel((int) v.getX(), (int) v.getY(), (int) v.getZ(), (float) nextLabel)
        }
    }

    splitRecords << [preLabel  : preLabel,
                     oldId     : pid,
                     parentVox : rootVox,
                     ratio     : rootVox / preMedianVolume,
                     nLevels   : levels.size(),
                     nParts    : leaves.size(),
                     damVoxels : totalDam,
                     parts     : partInfos,
                     levels    : levels,
                     bbox      : [x0: object.getXmin(), x1: object.getXmax(),
                                  y0: object.getYmin(), y1: object.getYmax()]]

    log("  split candidate v1.0 #${pid}: ${rootVox} vox " +
        "(${String.format('%.2f', rootVox / preMedianVolume)} x median) -> " +
        "${leaves.size()} final part(s) after ${levels.size() - 1} round(s), " +
        "${totalDam} dam voxels lost total" + (leaves.size() <= 1 ? " [NOT SPLIT]" : ""))
}
log("watershed split done: ${splitRecords.count { it.nParts > 1 }} object(s) actually split, " +
    "${splitRecords.count { it.nParts <= 1 }} left intact; " +
    "${splitEvents.size()} split event(s) over ${splitEvents.collect { it.depth }.unique().sort()} depth level(s); " +
    "total labels now ${nextLabel}")

// --- 4) 拆分后重新建对象、重新做体积筛选 ---
Objects3DPopulation population = new Objects3DPopulation(nucleusLabels, 0)
population.setCalibration(sx, sz, cal.getUnit())
def allObjects = population.getObjectsList()
def kept = allObjects.findAll { it.getVolumePixels() >= MIN_NUCLEUS_VOXELS }
int rejected = allObjects.size() - kept.size()
log("post-split objects: ${allObjects.size()}; passing volume (>= ${MIN_NUCLEUS_VOXELS} voxels): " +
    "${kept.size()}; rejected: ${rejected}")
if (kept.isEmpty()) {
    imp.close()
    throw new IllegalStateException("No nucleus passed MIN_NUCLEUS_VOXELS=${MIN_NUCLEUS_VOXELS}")
}

// 两张表：
//   <stem>_split_objects.csv  每一次真正发生的拆分一组行，带 depth，看得出刀是第几轮下的
//   <stem>_split_leaves.csv   递归结束后每个母对象的最终块，带该块是第几层产生的
def splitCsv = new File(REPORT_DIR, stem + "_split_objects.csv")
splitCsv.text = "old_cell_id_v1_0,depth,parent_part_id,parent_part_volume_voxels," +
        "parent_over_pre_median,root_volume_voxels,n_seeds,n_parts," +
        "part_id,part_volume_voxels,part_fraction_of_parent,part_fraction_of_root," +
        "part_still_over_cutoff,dam_voxels_lost\n"
splitEvents.each { ev ->
    ev.parts.each { p ->
        splitCsv << [ev.oldId, ev.depth, ev.parentId, ev.parentVox,
                     String.format("%.3f", ev.parentVox / preMedianVolume), ev.rootVox,
                     ev.nSeeds, ev.nParts,
                     p.id, p.voxels,
                     String.format("%.4f", p.voxels / (double) ev.parentVox),
                     String.format("%.4f", p.voxels / (double) ev.rootVox),
                     p.voxels > splitVolumeCutoff ? 1 : 0,
                     ev.dam].join(",") + "\n"
    }
}
def leavesCsv = new File(REPORT_DIR, stem + "_split_leaves.csv")
leavesCsv.text = "old_cell_id_v1_0,root_volume_voxels,root_over_pre_median,n_rounds,n_final_parts," +
        "part_id,part_depth,part_volume_voxels,part_over_pre_median,part_fraction_of_root," +
        "part_passes_min_voxels,total_dam_voxels_lost\n"
splitRecords.each { rec ->
    rec.parts.each { p ->
        leavesCsv << [rec.oldId, rec.parentVox, String.format("%.3f", rec.ratio),
                      rec.nLevels - 1, rec.nParts,
                      p.oldIdText, p.depth, p.voxels,
                      String.format("%.3f", p.voxels / preMedianVolume),
                      String.format("%.4f", p.voxels / (double) rec.parentVox),
                      p.voxels >= MIN_NUCLEUS_VOXELS ? 1 : 0,
                      rec.damVoxels].join(",") + "\n"
    }
}
log("split tables written to ${splitCsv.name} (${splitEvents.size()} split event(s)) and " +
    "${leavesCsv.name} (${splitRecords.size()} candidate object(s), " +
    "${splitRecords.sum(0) { it.nParts }} final part(s))")

def volumes = kept.collect { (double) it.getVolumePixels() }.sort()
double medianVolume = medianOf(volumes)
log("post-split nucleus volume voxels: min=${(int) volumes.first()} median=${(int) medianVolume} " +
    "max=${(int) volumes.last()}")

// 显示/遍历顺序：在拆分后的全集上按行优先排一次。
def orderedAll = rowMajorOrder(kept, ROW_BAND_PX)
def oldIdOf = [:]     // new label -> 显示用 v1.0 编号文本
kept.each { oldIdOf[it.getValue()] = oldIdTextOf[it.getValue()] }

// 每个核的裁切框。中心体的局部阈值窗口、实际输出的 tif、污染检查，全部用这同一个框。
def cropBoxOf = [:]
kept.each { Object3D object ->
    int x0 = Math.max(0, object.getXmin() - CROP_MARGIN_PX)
    int y0 = Math.max(0, object.getYmin() - CROP_MARGIN_PX)
    int x1 = Math.min(width - 1, object.getXmax() + CROP_MARGIN_PX)
    int y1 = Math.min(height - 1, object.getYmax() + CROP_MARGIN_PX)
    cropBoxOf[object.getValue()] = [x0: x0, y0: y0, x1: x1, y1: y1]
}

// --- 5) 中心体候选：逐核局部 MaxEntropy（沿用 1.2）---
// 亮度一律从未滤波的 C3 读，与基线一致。
ImagePlus spotDisplay = extractChannel(imp, SPOT_CHANNEL)
ImageHandler spotIntensity = ImageHandler.wrap(spotDisplay)

// 5a) 1.1 的全局阈值：只算出来做对照分类，不参与本版的检出。
ImagePlus spotsGlobalBlur = extractChannel(imp, SPOT_CHANNEL)
GaussianBlur3D.blur(spotsGlobalBlur, SPOT_SIGMA_XY, SPOT_SIGMA_XY, SPOT_SIGMA_Z)
double globalSpotThreshold = autoThreshold(spotsGlobalBlur, SPOT_THRESHOLD_METHOD)
spotsGlobalBlur.close()
log("reference-only: v1.1 whole-image ${SPOT_THRESHOLD_METHOD} threshold = ${globalSpotThreshold}")

// 5b) 每个核的窗口各算一次 MaxEntropy，二值结果 OR 进同一张全图掩膜。
ImageByte spotMaskGlobal = new ImageByte("spot mask", width, height, slices)
def localThresholdOf = [:]
orderedAll.each { Object3D nucleus ->
    def box = cropBoxOf[nucleus.getValue()]
    int w = box.x1 - box.x0 + 1
    int h = box.y1 - box.y0 + 1
    ImagePlus sub = cropChannel(imp, SPOT_CHANNEL, box.x0, box.y0, w, h)
    GaussianBlur3D.blur(sub, SPOT_SIGMA_XY, SPOT_SIGMA_XY, SPOT_SIGMA_Z)
    double t = autoThreshold(sub, SPOT_THRESHOLD_METHOD)
    localThresholdOf[nucleus.getValue()] = t
    ImageStack subStack = sub.getStack()
    for (int z = 0; z < slices; z++) {
        ImageProcessor ip = subStack.getProcessor(z + 1)
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (ip.getf(x, y) >= t) spotMaskGlobal.setPixel(box.x0 + x, box.y0 + y, z, 255.0f)
            }
        }
    }
    sub.close()
}
def thresholdValues = localThresholdOf.values().sort()
log("per-nucleus local ${SPOT_THRESHOLD_METHOD} thresholds over ${thresholdValues.size()} windows: " +
    "min=${String.format('%.1f', thresholdValues.first())} " +
    "median=${String.format('%.1f', thresholdValues[(int) (thresholdValues.size() / 2)])} " +
    "max=${String.format('%.1f', thresholdValues.last())} " +
    "(v1.1 global was ${String.format('%.1f', globalSpotThreshold)})")

ImageLabeller spotLabeller = new ImageLabeller()
spotLabeller.setMinSize(MIN_SPOT_VOXELS)
ImageInt spotLabels = spotLabeller.getLabels(spotMaskGlobal, false)
Objects3DPopulation spotPopulation = new Objects3DPopulation(spotLabels, 0)
spotPopulation.setCalibration(sx, sz, cal.getUnit())
def allSpots = spotPopulation.getObjectsList()
def spotCandidates = allSpots.findAll {
    it.getVolumePixels() >= MIN_SPOT_VOXELS && it.getVolumePixels() <= MAX_SPOT_VOXELS
}
int oversizedSpots = allSpots.count { it.getVolumePixels() > MAX_SPOT_VOXELS }
log("spot objects >= ${MIN_SPOT_VOXELS} voxels: ${allSpots.size()}; " +
    "candidates in [${MIN_SPOT_VOXELS}, ${MAX_SPOT_VOXELS}]: ${spotCandidates.size()}; " +
    "oversized rejected: ${oversizedSpots}")

// 5c) 对照报告：有多少候选的 maxI 低于 1.1 的全局阈值（= 1.1 漏掉的那些）。
def belowGlobal = spotCandidates.findAll { it.getPixMaxValue(spotIntensity) < globalSpotThreshold }
log("candidates whose max intensity is BELOW the v1.1 global threshold " +
    "(${String.format('%.1f', globalSpotThreshold)}): ${belowGlobal.size()} of ${spotCandidates.size()}")
int belowGate = spotCandidates.count { it.getPixMaxValue(spotIntensity) < MIN_SPOT_MAX_INTENSITY }
log("candidates whose max intensity is BELOW the v1.3 intensity gate " +
    "(${MIN_SPOT_MAX_INTENSITY}): ${belowGate} of ${spotCandidates.size()}")

// --- 6) 归属判定：每个候选分给“离它最近的那个核”（质心到核表面的标定距离）---
// 精确最近邻，带外接框下界剪枝；结果与全算一遍逐字相同。
def assignedOf = [:]
def distanceOf = [:]
def candidatesOfNucleus = [:]
kept.each { candidatesOfNucleus[it.getValue()] = [] }

spotCandidates.each { Object3D spot ->
    double cx = spot.getCenterX(), cy = spot.getCenterY(), cz = spot.getCenterZ()
    def ranked = kept.collect { Object3D nucleus ->
        [nucleus: nucleus, bound: bboxLowerBoundUnit(nucleus, cx, cy, cz, sx, sy, sz)]
    }.sort { it.bound }
    Object3D best = null
    double bestDistance = Double.MAX_VALUE
    for (entry in ranked) {
        if (entry.bound >= bestDistance) break
        double d = spot.distCenterBorderUnit(entry.nucleus)
        if (d < bestDistance) { bestDistance = d; best = entry.nucleus }
    }
    if (best != null) {
        assignedOf[spot.getValue()] = best
        distanceOf[spot.getValue()] = bestDistance
        candidatesOfNucleus[best.getValue()] << spot
    }
}

def spotsCsv = new File(REPORT_DIR, stem + "_spot_candidates.csv")
spotsCsv.text = "spot_label,volume_voxels,volume_um3,max_intensity,mean_intensity," +
        "centroid_x_px,centroid_y_px,centroid_z_px," +
        "assigned_nucleus_label,assigned_old_cell_id,distance_to_that_nucleus_surface_um," +
        "within_${(int) MAX_OWN_SPOT_DISTANCE_UM}um,passes_intensity_gate," +
        "centroid_inside_some_nucleus_label,below_v11_global_threshold\n"
spotCandidates.sort(false) { -it.getPixMaxValue(spotIntensity) }.each { Object3D spot ->
    Object3D nucleus = assignedOf[spot.getValue()]
    double d = distanceOf[spot.getValue()] == null ? Double.NaN : distanceOf[spot.getValue()]
    int insideLabel = labelAt(nucleusLabels, (int) Math.round(spot.getCenterX()),
            (int) Math.round(spot.getCenterY()), (int) Math.round(spot.getCenterZ()), width, height, slices)
    spotsCsv << [spot.getValue(), (int) spot.getVolumePixels(),
                 String.format("%.3f", spot.getVolumePixels() * voxelVolume),
                 String.format("%.0f", spot.getPixMaxValue(spotIntensity)),
                 String.format("%.1f", spot.getPixMeanValue(spotIntensity)),
                 String.format("%.2f", spot.getCenterX()), String.format("%.2f", spot.getCenterY()),
                 String.format("%.2f", spot.getCenterZ()),
                 nucleus == null ? "" : nucleus.getValue(),
                 nucleus == null ? "" : oldIdOf[nucleus.getValue()],
                 Double.isNaN(d) ? "" : String.format("%.4f", d),
                 (!Double.isNaN(d) && d <= MAX_OWN_SPOT_DISTANCE_UM) ? 1 : 0,
                 spot.getPixMaxValue(spotIntensity) >= MIN_SPOT_MAX_INTENSITY ? 1 : 0,
                 insideLabel,
                 spot.getPixMaxValue(spotIntensity) < globalSpotThreshold ? 1 : 0].join(",") + "\n"
}
log("spot candidate table written to ${spotsCsv.name}")

// --- 7) 保留条件（1.3 的改动 1）---
// 取该核分到的候选中最亮的那个；它必须 maxI >= MIN_SPOT_MAX_INTENSITY 且距离 <= 5 um。
def ownSpotOf = [:]
def ownDistanceOf = [:]
def nearestSpotOf = [:]
def keptCells = []
def droppedCells = []
def dropReasonOf = [:]
def manualOverrideOf = [:]
def manualKeepOldIds = stem == "4" ? MANUAL_KEEP_SERIES_4 : (stem == "5" ? MANUAL_KEEP_SERIES_5 : [] as Set)
double postMedianVolume = medianOf(kept.collect { (double) it.getVolumePixels() })

// 按实际输出裁切框预先检查邻核和邻近中心体。先用 XY 外接框筛掉绝大多数不相交对象，
// 再仅遍历可能相交核的体素；到达阈值就停止计数。这样保留轻微边缘碎片，剔除真实邻核。
def cropContaminationOf = [:]
kept.each { Object3D nucleus ->
    def box = cropBoxOf[nucleus.getValue()]
    def foreignSpots = spotCandidates.findAll { Object3D spot ->
        Object3D owner = assignedOf[spot.getValue()]
        owner != null && owner.getValue() != nucleus.getValue() &&
                inRect(spot.getCenterX(), spot.getCenterY(), box.x0, box.y0, box.x1, box.y1)
    }
    def foreignNucleusCounts = [:]
    kept.each { Object3D other ->
        if (other.getValue() == nucleus.getValue() ||
                !bboxIntersectsRect(other, box.x0, box.y0, box.x1, box.y1)) return
        // 记录完整实际重叠量；不能只记到 1000，否则 CSV 会丢失清理体素数。
        int overlap = countVoxelsInRectUpTo(other, box.x0, box.y0, box.x1, box.y1,
                                            Integer.MAX_VALUE)
        if (overlap >= MIN_FOREIGN_NUCLEUS_VOXELS_IN_CROP) {
            foreignNucleusCounts[other.getValue()] = overlap
        }
    }
    def flags = []
    if (!foreignNucleusCounts.isEmpty()) flags << "裁切框含邻近细胞核 >= ${MIN_FOREIGN_NUCLEUS_VOXELS_IN_CROP} 体素"
    if (!foreignSpots.isEmpty() && !foreignNucleusCounts.isEmpty()) flags << "裁切框含邻近中心体候选"
    cropContaminationOf[nucleus.getValue()] = [flags: flags, foreignNucleusCounts: foreignNucleusCounts]
}
kept.each { Object3D nucleus ->
    def mine = candidatesOfNucleus[nucleus.getValue()]
    Object3D brightest = mine.isEmpty() ? null : mine.max { it.getPixMaxValue(spotIntensity) }
    Object3D nearest = mine.isEmpty() ? null : mine.min { distanceOf[it.getValue()] }
    double d = brightest == null ? Double.NaN : distanceOf[brightest.getValue()]
    ownSpotOf[nucleus.getValue()] = brightest
    ownDistanceOf[nucleus.getValue()] = d
    nearestSpotOf[nucleus.getValue()] = nearest
    if (brightest == null) {
        droppedCells << nucleus
        dropReasonOf[nucleus.getValue()] = "未分配到本核的中心体候选"
        return
    }
    double maxI = brightest.getPixMaxValue(spotIntensity)
    boolean brightEnough = maxI >= MIN_SPOT_MAX_INTENSITY
    boolean closeEnough = d <= MAX_OWN_SPOT_DISTANCE_UM
    boolean unresolvedFused = nucleus.getVolumePixels() > FUSED_VOLUME_RATIO * postMedianVolume
    boolean manuallyKept = manualKeepOldIds.contains(oldIdOf[nucleus.getValue()])
    def why = []
    if (!brightEnough) why << "本核最亮中心体候选强度 ${(int) maxI} < ${(int) MIN_SPOT_MAX_INTENSITY}"
    if (!closeEnough) why << "本核最亮中心体候选距离 ${String.format('%.3f', d)} um > ${MAX_OWN_SPOT_DISTANCE_UM} um"
    if (unresolvedFused) why << "分水岭后核体积仍为中位数 ${String.format('%.2f', nucleus.getVolumePixels() / postMedianVolume)} 倍，疑似粘连核"
    if (REQUIRE_CLEAN_SINGLE_CELL_CROP) why.addAll(cropContaminationOf[nucleus.getValue()].flags)
    if (manuallyKept) {
        keptCells << nucleus
        manualOverrideOf[nucleus.getValue()] = "用户目检确认保留"
        log("  MANUAL KEEP v1.0 #${oldIdOf[nucleus.getValue()]}: user-confirmed despite automatic flags: ${why.join('；')}")
    } else if (why.isEmpty()) {
        keptCells << nucleus
    } else {
        droppedCells << nucleus
        dropReasonOf[nucleus.getValue()] = why.join("；")
    }
}
log("nuclei kept (automatic or user-confirmed manual keep): ${keptCells.size()}; manual keeps=${manualOverrideOf.values().size()}")
log("nuclei dropped (no own centrosome, contamination, or unresolved fusion): ${droppedCells.size()} " +
    "-> v1.0 ids ${droppedCells.collect { oldIdOf[it.getValue()] }.sort()}")

def nucleiCsv = new File(REPORT_DIR, stem + "_nuclei_all.csv")
nucleiCsv.text = "old_cell_id_v1_0,new_cell_id_v1_5_2,nucleus_label,nucleus_volume_voxels,was_split," +
        "centroid_x_px,centroid_y_px,centroid_z_px,local_spot_threshold,n_assigned_candidates," +
        "brightest_own_candidate_um,brightest_own_candidate_voxels,brightest_own_candidate_max_intensity," +
        "nearest_own_candidate_um,nearest_own_candidate_voxels,nearest_own_candidate_max_intensity," +
        "nearest_any_candidate_um,nearest_any_candidate_voxels,nearest_any_candidate_max_intensity," +
        "touches_xy_edge,kept,manual_override,drop_reason\n"

// --- 8) 剔除后重新行优先编号（规则与 1.0 相同，作用在保留子集上）---
def ordered = rowMajorOrder(keptCells, ROW_BAND_PX)
def newIdOf = [:]
ordered.eachWithIndex { Object3D object, int index -> newIdOf[object.getValue()] = index + 1 }
log("renumbered rows: ${ordered.size()} cells")

// --- 9) 逐个裁切全部通道 / 全部 z 层并存 tif ---
def cellsCsv = new File(REPORT_DIR, stem + "_cells.csv")
cellsCsv.text = "cell_id,old_cell_id_v1_0,file_name,crop_x0,crop_y0,crop_width_px,crop_height_px," +
        "nucleus_bbox_x0,nucleus_bbox_x1,nucleus_bbox_y0,nucleus_bbox_y1,nucleus_bbox_z0,nucleus_bbox_z1," +
        "nucleus_volume_voxels,nucleus_volume_um3,volume_over_median,was_split," +
        "centroid_x_px,centroid_y_px,centroid_z_px,local_spot_threshold," +
        "own_spot_distance_um,own_spot_voxels,own_spot_max_intensity,own_spot_below_v11_global,n_own_candidates," +
        "nearest_spot_distance_um,nearest_spot_max_intensity," +
        "own_spot_in_crop,foreign_spots_in_crop,foreign_spot_owner_old_ids," +
        "foreign_nucleus_voxels_in_crop,foreign_nucleus_old_ids,contamination_flags," +
        "foreign_nucleus_cleaned,foreign_spot_objects_cleaned," +
        "touches_xy_edge,margin_clipped,suspected_fused,manual_override,channels,slices,bit_depth," +
        "pixel_width_um,pixel_height_um,voxel_depth_um\n"

def splitLabels = new HashSet()
splitRecords.findAll { it.nParts > 1 }.each { rec -> rec.parts.each { splitLabels << it.label } }

def records = []
ordered.eachWithIndex { Object3D object, int index ->
    int id = index + 1
    String oldId = oldIdOf[object.getValue()]
    String fileName = namePrefix + id + ".tif"
    int nx0 = object.getXmin(); int nx1 = object.getXmax()
    int ny0 = object.getYmin(); int ny1 = object.getYmax()
    def box = cropBoxOf[object.getValue()]
    int x0 = box.x0, y0 = box.y0, x1 = box.x1, y1 = box.y1
    int cropW = x1 - x0 + 1
    int cropH = y1 - y0 + 1
    boolean touchesEdge = (nx0 <= 0 || ny0 <= 0 || nx1 >= width - 1 || ny1 >= height - 1)
    boolean marginClipped = (nx0 - CROP_MARGIN_PX < 0 || ny0 - CROP_MARGIN_PX < 0 ||
            nx1 + CROP_MARGIN_PX > width - 1 || ny1 + CROP_MARGIN_PX > height - 1)
    int volume = (int) object.getVolumePixels()
    boolean fused = volume > FUSED_VOLUME_RATIO * medianVolume
    boolean wasSplit = splitLabels.contains(object.getValue())
    boolean manualOverride = manualOverrideOf.containsKey(object.getValue())

    Object3D ownSpot = ownSpotOf[object.getValue()]
    double ownDistance = ownDistanceOf[object.getValue()]
    Object3D nearSpot = nearestSpotOf[object.getValue()]
    int nOwn = candidatesOfNucleus[object.getValue()].size()

    def foreignSpots = spotCandidates.findAll { Object3D spot ->
        Object3D owner = assignedOf[spot.getValue()]
        owner != null && owner.getValue() != object.getValue() &&
                inRect(spot.getCenterX(), spot.getCenterY(), x0, y0, x1, y1)
    }
    def foreignSpotOwners = foreignSpots.collect { oldIdOf[assignedOf[it.getValue()].getValue()] }
                                        .unique().sort()
    boolean ownSpotInCrop = ownSpot != null &&
            inRect(ownSpot.getCenterX(), ownSpot.getCenterY(), x0, y0, x1, y1)

    // 1.4.2 复用实际体素重叠预检结果；保留对象在此阈值下应无污染。
    def foreignNucleusCounts = cropContaminationOf[object.getValue()].foreignNucleusCounts
    int foreignNucleusVoxelCount = (int) foreignNucleusCounts.values().sum(0)
    def foreignNucleusIds = foreignNucleusCounts.keySet().collect { oldIdOf[it] }.sort()

    def flags = []
    if (!ownSpotInCrop) flags << "OWN_SPOT_OUTSIDE_CROP"
    if (!foreignSpots.isEmpty()) flags << "FOREIGN_SPOT"
    if (foreignNucleusVoxelCount >= MIN_FOREIGN_NUCLEUS_VOXELS_IN_CROP) flags << "FOREIGN_NUCLEUS"

    ImagePlus crop = cropAllChannelsCleanForeign(imp, x0, y0, cropW, cropH,
            nucleusLabels, object.getValue(), spotLabels, ownSpot.getValue(),
            NUCLEUS_CHANNEL, SPOT_CHANNEL, displayReference)
    if (SAVE_TIFS) {
        new FileSaver(crop).saveAsTiff(new File(OUTPUT_DIR, fileName).absolutePath)
    }

    records << [id: id, oldId: oldId, name: fileName, x0: x0, y0: y0, w: cropW, h: cropH,
                volume: volume, fused: fused, edge: touchesEdge, split: wasSplit, manual: manualOverride,
                distance: ownDistance, maxI: (int) ownSpot.getPixMaxValue(spotIntensity),
                flags: flags,
                blue: MAKE_OVERVIEW ? maxProjection(crop, NUCLEUS_CHANNEL) : null,
                red: MAKE_OVERVIEW ? maxProjection(crop, SPOT_CHANNEL) : null]

    cellsCsv << [id, oldId, fileName, x0, y0, cropW, cropH,
                 nx0, nx1, ny0, ny1, object.getZmin(), object.getZmax(),
                 volume, String.format("%.3f", volume * voxelVolume),
                 String.format("%.3f", volume / medianVolume), wasSplit ? 1 : 0,
                 String.format("%.2f", object.getCenterX()), String.format("%.2f", object.getCenterY()),
                 String.format("%.2f", object.getCenterZ()),
                 String.format("%.1f", localThresholdOf[object.getValue()]),
                 String.format("%.4f", ownDistance), (int) ownSpot.getVolumePixels(),
                 String.format("%.0f", ownSpot.getPixMaxValue(spotIntensity)),
                 ownSpot.getPixMaxValue(spotIntensity) < globalSpotThreshold ? 1 : 0, nOwn,
                 nearSpot == null ? "" : String.format("%.4f", distanceOf[nearSpot.getValue()]),
                 nearSpot == null ? "" : String.format("%.0f", nearSpot.getPixMaxValue(spotIntensity)),
                 ownSpotInCrop ? 1 : 0, foreignSpots.size(),
                 foreignSpotOwners.isEmpty() ? "" : foreignSpotOwners.join(" "),
                 foreignNucleusVoxelCount,
                 foreignNucleusIds.isEmpty() ? "" : foreignNucleusIds.join(" "),
                 flags.isEmpty() ? "" : flags.join(" "),
                 foreignNucleusVoxelCount > 0 ? 1 : 0, foreignSpots.size(),
                 touchesEdge ? 1 : 0, marginClipped ? 1 : 0, fused ? 1 : 0, manualOverride ? 1 : 0,
                 crop.getNChannels(), crop.getNSlices(), crop.getBitDepth(),
                 crop.getCalibration().pixelWidth, crop.getCalibration().pixelHeight,
                 crop.getCalibration().pixelDepth].join(",") + "\n"
    log("${fileName} (v1.0 ${namePrefix}${oldId}): crop ${cropW}x${cropH} at (${x0},${y0}), " +
        "nucleus ${volume} vox, centrosome ${String.format('%.3f', ownDistance)} um " +
        "(${(int) ownSpot.getVolumePixels()} vox, maxI ${(int) ownSpot.getPixMaxValue(spotIntensity)}" +
        (ownSpot.getPixMaxValue(spotIntensity) < globalSpotThreshold ? ", MISSED BY v1.1" : "") + "), " +
        "${nOwn} assigned" +
        (wasSplit ? " [WATERSHED-SPLIT]" : "") +
        (manualOverride ? " [MANUAL-KEEP]" : "") +
        (touchesEdge ? " [XY-EDGE]" : "") + (fused ? " [FUSED?]" : "") +
        (flags.isEmpty() ? "" : " [" + flags.join("][") + "]"))
    crop.close()
}

// --- 10) 所有核的判定表（保留 + 剔除）---
def droppedRecords = []
orderedAll.each { Object3D nucleus ->
    int label = nucleus.getValue()
    Object3D ownSpot = ownSpotOf[label]
    Object3D nearSpot = nearestSpotOf[label]
    double ownDistance = ownDistanceOf[label]
    int nOwn = candidatesOfNucleus[label].size()
    boolean isKept = newIdOf.containsKey(label)
    Object3D nearestAny = spotCandidates.isEmpty() ? null :
            spotCandidates.min { it.distCenterBorderUnit(nucleus) }
    double nearestAnyDistance = nearestAny == null ? Double.NaN : nearestAny.distCenterBorderUnit(nucleus)
    boolean touchesEdge = (nucleus.getXmin() <= 0 || nucleus.getYmin() <= 0 ||
            nucleus.getXmax() >= width - 1 || nucleus.getYmax() >= height - 1)
    String reason = isKept ? "" : dropReasonOf[label]
    nucleiCsv << [oldIdOf[label], isKept ? newIdOf[label] : "", label, (int) nucleus.getVolumePixels(),
                  splitLabels.contains(label) ? 1 : 0,
                  String.format("%.2f", nucleus.getCenterX()), String.format("%.2f", nucleus.getCenterY()),
                  String.format("%.2f", nucleus.getCenterZ()),
                  String.format("%.1f", localThresholdOf[label]), nOwn,
                  Double.isNaN(ownDistance) ? "" : String.format("%.4f", ownDistance),
                  ownSpot == null ? "" : (int) ownSpot.getVolumePixels(),
                  ownSpot == null ? "" : String.format("%.0f", ownSpot.getPixMaxValue(spotIntensity)),
                  nearSpot == null ? "" : String.format("%.4f", distanceOf[nearSpot.getValue()]),
                  nearSpot == null ? "" : (int) nearSpot.getVolumePixels(),
                  nearSpot == null ? "" : String.format("%.0f", nearSpot.getPixMaxValue(spotIntensity)),
                  Double.isNaN(nearestAnyDistance) ? "" : String.format("%.4f", nearestAnyDistance),
                  nearestAny == null ? "" : (int) nearestAny.getVolumePixels(),
                  nearestAny == null ? "" : String.format("%.0f", nearestAny.getPixMaxValue(spotIntensity)),
                  touchesEdge ? 1 : 0, isKept ? 1 : 0,
                  manualOverrideOf.containsKey(label) ? manualOverrideOf[label] : "", reason].join(",") + "\n"

    if (!isKept && MAKE_DROPPED_PREVIEW) {
        int ww = Math.min(DROPPED_WINDOW_PX, width)
        int hh = Math.min(DROPPED_WINDOW_PX, height)
        int wx = (int) Math.round(nucleus.getCenterX()) - (int) (ww / 2)
        int wy = (int) Math.round(nucleus.getCenterY()) - (int) (hh / 2)
        wx = Math.max(0, Math.min(width - ww, wx))
        wy = Math.max(0, Math.min(height - hh, wy))
        // 窗口内的全部候选。1.4 只给其中一个画圆环（ringed=true）：优先是该核自己分到的
        // 最亮候选；该核一个都没分到时，改成窗口内最亮的那个，并标出它归属给谁。
        // 其余候选在图上只剩 2 px 小点，不带文字——1.3 给每个候选都画圆环加标签，
        // #19 那格 22 个圆环把细胞盖死了。
        def inWindow = spotCandidates.findAll {
            inRect(it.getCenterX(), it.getCenterY(), wx, wy, wx + ww - 1, wy + hh - 1)
        }
        Object3D ringSpot = ownSpot
        if (ringSpot == null || !inRect(ringSpot.getCenterX(), ringSpot.getCenterY(),
                                        wx, wy, wx + ww - 1, wy + hh - 1)) {
            ringSpot = inWindow.isEmpty() ? null : inWindow.max { it.getPixMaxValue(spotIntensity) }
        }
        def marks = inWindow.collect { Object3D spot ->
            Object3D owner = assignedOf[spot.getValue()]
            [x     : (int) Math.round(spot.getCenterX()) - wx,
             y     : (int) Math.round(spot.getCenterY()) - wy,
             owner : owner == null ? "?" : ("" + oldIdOf[owner.getValue()]),
             mine  : owner != null && owner.getValue() == label,
             bright: spot.getPixMaxValue(spotIntensity) >= MIN_SPOT_MAX_INTENSITY,
             ringed: ringSpot != null && spot.getValue() == ringSpot.getValue()]
        }
        droppedRecords << [oldId  : oldIdOf[label],
                           ordKey : (parentOldIdOf[label] == null ? 9999 : parentOldIdOf[label]),
                           volUm3 : nucleus.getVolumePixels() * voxelVolume,
                           near   : nearestAnyDistance,
                           maxI   : ownSpot == null ? -1 : (int) ownSpot.getPixMaxValue(spotIntensity),
                           reason : reason,
                           edge   : touchesEdge,
                           w      : ww, h: hh,
                           blue   : windowMaxProjection(imp, NUCLEUS_CHANNEL, wx, wy, ww, hh),
                           red    : windowMaxProjection(imp, SPOT_CHANNEL, wx, wy, ww, hh),
                           outline: borderOf(insideProjection(nucleusLabels, label, wx, wy, ww, hh, slices), ww, hh),
                           marks  : marks]
    }
}
log("per-nucleus decision table written to ${nucleiCsv.name}")

// --- 11) 总览图、被剔除细胞缩略图、拆分对照图 ---
if (MAKE_OVERVIEW) {
    File overview = new File(OUTPUT_DIR, stem + "_overview.png")
    saveOverview(records, overview, OVERVIEW_COLUMNS, OVERVIEW_PAD, OVERVIEW_HEADER,
                 OVERVIEW_MIN_TILE_W, OVERVIEW_CLIP_LOW, OVERVIEW_CLIP_HIGH)
    log("overview: " + overview.name + " (${records.size()} tiles) -> " + OUTPUT_DIR.name)
}
if (MAKE_DROPPED_PREVIEW && !droppedRecords.isEmpty()) {
    droppedRecords.sort { it.ordKey }
    File droppedPng = new File(DROPPED_DIR, stem + "_dropped.png")
    saveDroppedPreview(droppedRecords, droppedPng, stem, DROPPED_COLUMNS, OVERVIEW_PAD,
                       DROPPED_MARKER_RADIUS_PX, DROPPED_DOT_RADIUS_PX, DROPPED_TILE_MIN_W,
                       OVERVIEW_CLIP_LOW, OVERVIEW_CLIP_HIGH, (int) MIN_SPOT_MAX_INTENSITY)
    log("dropped preview: " + droppedPng.name + " (${droppedRecords.size()} tiles) -> " + DROPPED_DIR.name)
}
if (MAKE_SPLIT_QC && !splitRecords.isEmpty()) {
    // 每个母对象一行，一行里每一轮一格：第 0 格是拆分前，第 k 格是第 k 轮之后。
    // 每一块的颜色按它第一次出现的顺序分配，块在后续格里保持同色，所以能一眼看出
    // 哪一块在第几轮被切开、切成了哪两块。
    def qcTiles = splitRecords.collect { rec ->
        int wx = Math.max(0, rec.bbox.x0 - SPLIT_QC_MARGIN_PX)
        int wy = Math.max(0, rec.bbox.y0 - SPLIT_QC_MARGIN_PX)
        int wx1 = Math.min(width - 1, rec.bbox.x1 + SPLIT_QC_MARGIN_PX)
        int wy1 = Math.min(height - 1, rec.bbox.y1 + SPLIT_QC_MARGIN_PX)
        int ww = wx1 - wx + 1, hh = wy1 - wy + 1
        def colorIndexOf = [:]
        def panels = rec.levels.collect { level ->
            level.collect { piece ->
                if (!colorIndexOf.containsKey(piece.id)) colorIndexOf[piece.id] = colorIndexOf.size()
                [id     : piece.id,
                 voxels : piece.voxels.size(),
                 depth  : piece.depth,
                 colorIx: colorIndexOf[piece.id],
                 inside : projectionOfVoxels(piece.voxels, wx, wy, ww, hh)]
            }
        }
        [rec: rec, w: ww, h: hh,
         blue  : windowMaxProjection(imp, NUCLEUS_CHANNEL, wx, wy, ww, hh),
         red   : windowMaxProjection(imp, SPOT_CHANNEL, wx, wy, ww, hh),
         panels: panels]
    }
    File qcPng = new File(QC_DIR, stem + "_split_qc.png")
    saveSplitQc(qcTiles, qcPng, stem, OVERVIEW_PAD, OVERVIEW_CLIP_LOW, OVERVIEW_CLIP_HIGH,
                (int) preMedianVolume, SPLIT_VOLUME_RATIO, SPLIT_NOISE_TOLERANCE_UM)
    log("split qc: " + qcPng.name + " (${qcTiles.size()} candidate object(s), " +
        "${qcTiles.collect { it.panels.size() }.max()} panel(s) per row) -> " + QC_DIR.name)
}
if (MAKE_DECISION_MAP) {
    File decisionPng = new File(QC_DIR, stem + "_decision_map.png")
    saveDecisionMap(imp, nucleusLabels, orderedAll, newIdOf, decisionPng, stem,
                    OVERVIEW_PAD, OVERVIEW_CLIP_LOW, OVERVIEW_CLIP_HIGH,
                    NUCLEUS_CHANNEL, SPOT_CHANNEL)
    log("original-field decision map: " + decisionPng.name + " (green kept=${keptCells.size()}, yellow dropped=${droppedCells.size()}) -> " + QC_DIR.name)
}
log("figure text: ${TRUNCATED_STRINGS} string(s) were wider than their cell and got truncated " +
    "with \"...\" (0 means every label fits as drawn)")

spotIntensity.closeImagePlus()
nucleiRaw.close()
nucleiMask.close()
imp.close()

int edgeCount = records.count { it.edge }
int fusedCount = records.count { it.fused }
int flaggedCount = records.count { !it.flags.isEmpty() }
int splitKeptCount = records.count { it.split }
log("DONE: ${records.size()} all-channel/all-Z crops with selective foreign C1/C3 cleanup written (from ${kept.size()} volume-passing nuclei, " +
    "${splitRecords.count { it.nParts > 1 }} object(s) watershed-split), " +
    "${droppedCells.size()} dropped for no own centrosome, contamination, or unresolved fusion, " +
    "${edgeCount} touching XY edge, ${fusedCount} suspected fused, " +
    "${splitKeptCount} of the kept cells came from a watershed split, " +
    "${flaggedCount} with contamination flags, " +
    "${rejected} post-split components rejected by volume")

}   // end INPUT_FILES.each

displayReference.close()

// ===================== 辅助函数 =====================

// 递归拆分用的后缀字符：每一层换一套，看编号就知道这块是第几轮切出来的。
// 第 1 层 a b c...，第 2 层 1 2 3...，第 3 层 x y z...
String suffixChar(int depth, int k) {
    def sets = ["abcdefghij", "123456789", "xyzuvwrst"]
    String s = sets[(depth - 1) % sets.size()]
    return "" + s.charAt(Math.min(k, s.length() - 1))
}

// 对一组体素跑一次 EDT + MaximaFinder + Watershed3D。
// 返回 [nSeeds: 种子数, dam: 被分水岭吃掉的体素数, parts: 按体积降序的体素组列表]。
// parts.size() == 1 表示没拆开（这时 dam = 0，体素原样退回）。
Map watershedVoxels(List voxels, int pad, int W, int H, int D,
                    float sx, float sz, float radXY, float radZ, float tol, double wsThresh) {
    int x0 = Integer.MAX_VALUE, y0 = Integer.MAX_VALUE, z0 = Integer.MAX_VALUE
    int x1 = -1, y1 = -1, z1 = -1
    voxels.each { Voxel3D v ->
        int x = (int) v.getX(), y = (int) v.getY(), z = (int) v.getZ()
        if (x < x0) x0 = x; if (x > x1) x1 = x
        if (y < y0) y0 = y; if (y > y1) y1 = y
        if (z < z0) z0 = z; if (z > z1) z1 = z
    }
    x0 = Math.max(0, x0 - pad); y0 = Math.max(0, y0 - pad); z0 = Math.max(0, z0 - pad)
    x1 = Math.min(W - 1, x1 + pad); y1 = Math.min(H - 1, y1 + pad); z1 = Math.min(D - 1, z1 + pad)
    int sw = x1 - x0 + 1, sh = y1 - y0 + 1, sd = z1 - z0 + 1

    ImageByte subMask = new ImageByte("sub", sw, sh, sd)
    voxels.each { Voxel3D v ->
        subMask.setPixel((int) v.getX() - x0, (int) v.getY() - y0, (int) v.getZ() - z0, 255.0f)
    }
    // EDT：标定过的 3D 距离变换，值 = 该体素到块外最近点的距离（um），即内切球半径。
    def edt = EDT.run(subMask, 128.0f, sx, sz, false, 0)
    MaximaFinder maxima = new MaximaFinder(edt, radXY, radZ, tol)
    maxima.setVerbose(false)
    ImageHandler seeds = maxima.getImagePeaks()
    int nSeeds = maxima.getListPeaks() == null ? -1 : maxima.getListPeaks().size()
    if (nSeeds == 1) {
        // 只有一个极大，分水岭必然只还回一块，省掉一次泛洪。
        return [nSeeds: nSeeds, dam: 0, parts: [voxels]]
    }
    Watershed3D watershed = new Watershed3D(edt, seeds, wsThresh, 0)
    watershed.setLabelSeeds(true)
    ImageHandler wsImage = watershed.getWatershedImage3D()

    def byLabel = [:]
    int dam = 0
    voxels.each { Voxel3D v ->
        int wl = (int) wsImage.getPixel((int) v.getX() - x0, (int) v.getY() - y0, (int) v.getZ() - z0)
        if (wl == 0) { dam++; return }
        def list = byLabel[wl]
        if (list == null) { list = []; byLabel[wl] = list }
        list << v
    }
    def parts = byLabel.values().sort { -it.size() }
    if (parts.size() <= 1) return [nSeeds: nSeeds, dam: 0, parts: [voxels]]
    return [nSeeds: nSeeds, dam: dam, parts: parts]
}

// 一组体素在窗口内的 XY 投影。不需要先落进标签图，所以中间层的块也画得出来。
boolean[] projectionOfVoxels(List voxels, int x0, int y0, int w, int h) {
    boolean[] inside = new boolean[w * h]
    voxels.each { Voxel3D v ->
        int x = (int) v.getX() - x0
        int y = (int) v.getY() - y0
        if (x >= 0 && y >= 0 && x < w && y < h) inside[y * w + x] = true
    }
    return inside
}

double medianOf(List values) {
    def v = values.sort(false)
    return v.size() % 2 == 1 ? (double) v[(int) (v.size() / 2)] :
            0.5d * ((double) v[v.size().intdiv(2) - 1] + (double) v[v.size().intdiv(2)])
}

// 行优先排序：先按质心 y 分行（|y - 行锚点 y| <= band 算同一行），行内按 x 升序。
// 与 1.0 逐字相同。
List rowMajorOrder(List objects, double band) {
    def byY = objects.sort(false) { it.getCenterY() }
    def rows = []
    def currentRow = []
    double anchorY = Double.NaN
    byY.each { Object3D object ->
        if (currentRow.isEmpty() || object.getCenterY() - anchorY > band) {
            if (!currentRow.isEmpty()) rows << currentRow
            currentRow = []
            anchorY = object.getCenterY()
        }
        currentRow << object
    }
    if (!currentRow.isEmpty()) rows << currentRow
    def ordered = []
    rows.each { row -> ordered.addAll(row.sort(false) { it.getCenterX() }) }
    return ordered
}

// 点到一个对象外接框的标定距离：点到该对象表面真实距离的下界。
double bboxLowerBoundUnit(Object3D object, double x, double y, double z,
                          double sx, double sy, double sz) {
    double dx = Math.max(0.0d, Math.max(object.getXmin() - x, x - object.getXmax())) * sx
    double dy = Math.max(0.0d, Math.max(object.getYmin() - y, y - object.getYmax())) * sy
    double dz = Math.max(0.0d, Math.max(object.getZmin() - z, z - object.getZmax())) * sz
    return Math.sqrt(dx * dx + dy * dy + dz * dz)
}

boolean inRect(double x, double y, int x0, int y0, int x1, int y1) {
    return x >= x0 && x <= x1 && y >= y0 && y <= y1
}

// 对原始全 Z 裁切，XY 外接框相交就意味着裁切中可能保留第二个核的原始信号。
boolean bboxIntersectsRect(Object3D object, int x0, int y0, int x1, int y1) {
    return object.getXmax() >= x0 && object.getXmin() <= x1 &&
           object.getYmax() >= y0 && object.getYmin() <= y1
}

// 只数到判定阈值即可，避免每个裁切框重复扫描整个三维体积。
int countVoxelsInRectUpTo(Object3D object, int x0, int y0, int x1, int y1, int stopAt) {
    int count = 0
    for (Voxel3D voxel : object.getVoxels()) {
        double x = voxel.getX(), y = voxel.getY()
        if (x >= x0 && x <= x1 && y >= y0 && y <= y1) {
            count++
            if (count >= stopAt) return count
        }
    }
    return count
}

Set keptLabels(List nuclei) {
    def set = new HashSet()
    nuclei.each { set << it.getValue() }
    return set
}

int labelAt(ImageInt labels, int x, int y, int z, int w, int h, int d) {
    if (x < 0 || y < 0 || z < 0 || x >= w || y >= h || z >= d) return 0
    return (int) labels.getPixel(x, y, z)
}

Map foreignNucleusVoxels(ImageInt labels, int x0, int y0, int x1, int y1, int slices,
                         int ownLabel, Set validLabels) {
    def counts = [:]
    for (int z = 0; z < slices; z++) {
        for (int y = y0; y <= y1; y++) {
            for (int x = x0; x <= x1; x++) {
                int label = (int) labels.getPixel(x, y, z)
                if (label != 0 && label != ownLabel && validLabels.contains(label)) {
                    counts[label] = (counts[label] ?: 0) + 1
                }
            }
        }
    }
    return counts
}

// 某个标签在窗口内的 XY 投影（全 z 的 OR）。
boolean[] insideProjection(ImageInt labels, int label, int x0, int y0, int w, int h, int slices) {
    boolean[] inside = new boolean[w * h]
    for (int z = 0; z < slices; z++) {
        for (int y = 0; y < h; y++) {
            int base = y * w
            for (int x = 0; x < w; x++) {
                if (!inside[base + x] && (int) labels.getPixel(x0 + x, y0 + y, z) == label) {
                    inside[base + x] = true
                }
            }
        }
    }
    return inside
}

// 投影的边界像素：属于该投影、且 4 邻域中有一个不属于它。只标边界，不填充。
boolean[] borderOf(boolean[] inside, int w, int h) {
    boolean[] border = new boolean[w * h]
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int i = y * w + x
            if (!inside[i]) continue
            boolean edge = (x == 0 || y == 0 || x == w - 1 || y == h - 1) ||
                    !inside[i - 1] || !inside[i + 1] || !inside[i - w] || !inside[i + w]
            border[i] = edge
        }
    }
    return border
}

// 单通道全 z 栈，用于分割。与冻结基线 1.0 的同名函数一致。
ImagePlus extractChannel(ImagePlus source, int channel) {
    ImageStack sourceStack = source.getStack()
    int slices = source.getNSlices()
    int frames = source.getNFrames()
    ImageStack result = new ImageStack(source.getWidth(), source.getHeight())
    for (int t = 1; t <= frames; t++) {
        for (int z = 1; z <= slices; z++) {
            result.addSlice(sourceStack.getProcessor(source.getStackIndex(channel, z, t)).duplicate())
        }
    }
    ImagePlus output = new ImagePlus("C" + channel, result)
    output.setDimensions(1, slices, frames)
    output.setCalibration(source.getCalibration().copy())
    return output
}

// 单通道、XY 裁切、全部 z 层。与基线在单细胞 tif 上看到的数据完全相同：
// 基线是先裁后模糊后阈值，这里也是先裁后模糊后阈值。
ImagePlus cropChannel(ImagePlus source, int channel, int x0, int y0, int w, int h) {
    ImageStack sourceStack = source.getStack()
    int slices = source.getNSlices()
    ImageStack result = new ImageStack(w, h)
    for (int z = 1; z <= slices; z++) {
        ImageProcessor ip = sourceStack.getProcessor(source.getStackIndex(channel, z, 1))
        ip.setRoi(x0, y0, w, h)
        result.addSlice(ip.crop())
        ip.resetRoi()
    }
    ImagePlus output = new ImagePlus("C" + channel + " crop", result)
    output.setDimensions(1, slices, 1)
    output.setCalibration(source.getCalibration().copy())
    return output
}

double autoThreshold(ImagePlus source, AutoThresholder.Method method) {
    StackStatistics statistics = new StackStatistics(source)
    int bin = new AutoThresholder().getThreshold(method, intHistogram(statistics.histogram))
    return statistics.histMin + (bin + 0.5) * statistics.binSize
}

ImagePlus thresholdStack(ImagePlus source, double threshold) {
    ImageStack result = new ImageStack(source.getWidth(), source.getHeight())
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
    output.setDimensions(1, source.getNSlices(), source.getNFrames())
    output.setCalibration(source.getCalibration().copy())
    return output
}

int[] intHistogram(long[] histogram) {
    int[] converted = new int[histogram.length]
    for (int index = 0; index < histogram.length; index++) {
        converted[index] = (int) Math.min(Integer.MAX_VALUE, histogram[index])
    }
    return converted
}

// 裁 XY，保留全部通道、全部 z 层、原位深与原标定。切片顺序按 ImageJ 的 XYCZT。
ImagePlus cropAllChannels(ImagePlus source, int x0, int y0, int w, int h) {
    ImageStack sourceStack = source.getStack()
    int channels = source.getNChannels()
    int slices = source.getNSlices()
    int frames = source.getNFrames()
    ImageStack result = new ImageStack(w, h)
    for (int t = 1; t <= frames; t++) {
        for (int z = 1; z <= slices; z++) {
            for (int c = 1; c <= channels; c++) {
                ImageProcessor ip = sourceStack.getProcessor(source.getStackIndex(c, z, t))
                ip.setRoi(x0, y0, w, h)
                result.addSlice(ip.crop())
                ip.resetRoi()
            }
        }
    }
    ImagePlus output = new ImagePlus("crop", result)
    output.setDimensions(channels, slices, frames)
    output.setOpenAsHyperStack(true)
    output.setCalibration(source.getCalibration().copy())
    // Preserve the source display contract. The pixel data are already raw crops;
    // this only copies Composite/LUT state so Fiji opens C1/C2/C3/C4 with the
    // same colours and mode as the original multi-channel TIFF.
    CompositeImage composite = new CompositeImage(output, CompositeImage.COMPOSITE)
    if (source instanceof CompositeImage) composite.copyLuts(source)
    composite.setMode(CompositeImage.COMPOSITE)
    String info = source.getInfoProperty()
    if (info != null) composite.setProperty("Info", info)
    return composite
}

// 原始矩形裁切的最小结构清理：背景和非目标通道完全保持原值；仅在 C1 清除其他核标签，
// 在 C3 清除其他中心体标签。目标核/中心体像素逐值不变，仍保留所有通道和全部 Z 层。
ImagePlus cropAllChannelsCleanForeign(ImagePlus source, int x0, int y0, int w, int h,
                                      ImageInt nucleusLabels, int nucleusLabel,
                                      ImageInt spotLabels, int spotLabel,
                                      int nucleusChannel, int spotChannel,
                                      ImagePlus displayReference) {
    ImageStack sourceStack = source.getStack()
    int channels = source.getNChannels()
    int slices = source.getNSlices()
    int frames = source.getNFrames()
    ImageStack result = new ImageStack(w, h)
    for (int t = 1; t <= frames; t++) {
        for (int z = 1; z <= slices; z++) {
            for (int c = 1; c <= channels; c++) {
                int sourceIndex = source.getStackIndex(c, z, t)
                ImageProcessor src = sourceStack.getProcessor(sourceIndex)
                src.setRoi(x0, y0, w, h)
                ImageProcessor out = src.crop()
                if (c == nucleusChannel || c == spotChannel) {
                    for (int y = 0; y < h; y++) {
                        for (int x = 0; x < w; x++) {
                            int gx = x0 + x, gy = y0 + y, gz = z - 1
                            if (c == nucleusChannel) {
                                int label = (int) nucleusLabels.getPixel(gx, gy, gz)
                                if (label != 0 && label != nucleusLabel) out.setf(x, y, 0.0f)
                            } else {
                                int label = (int) spotLabels.getPixel(gx, gy, gz)
                                if (label != 0 && label != spotLabel) out.setf(x, y, 0.0f)
                            }
                        }
                    }
                }
                result.addSlice(sourceStack.getSliceLabel(sourceIndex), out)
                src.resetRoi()
            }
        }
    }
    ImagePlus output = new ImagePlus("clean crop", result)
    output.setDimensions(channels, slices, frames)
    output.setOpenAsHyperStack(true)
    output.setCalibration(source.getCalibration().copy())
    CompositeImage composite = new CompositeImage(output, CompositeImage.COMPOSITE)
    if (displayReference instanceof CompositeImage) composite.copyLuts(displayReference)
    composite.setMode(CompositeImage.COMPOSITE)
    String info = source.getInfoProperty()
    if (info != null) composite.setProperty("Info", info)
    return composite
}

// 结构掩膜裁切：TIFF 外形仍是矩形（格式要求），但有效像素只保留目标核和目标中心体。
// C1 只保留目标核标签，C3 只保留归属目标核的中心体标签，其他通道保留两者的并集。
ImagePlus cropMaskedTarget(ImagePlus source, int x0, int y0, int w, int h,
                           ImageInt nucleusLabels, int nucleusLabel,
                           ImageInt spotLabels, int spotLabel,
                           int nucleusChannel, int spotChannel) {
    int channels = source.getNChannels()
    int slices = source.getNSlices()
    int frames = source.getNFrames()
    ImageStack result = new ImageStack(w, h)
    for (int t = 1; t <= frames; t++) {
        for (int z = 1; z <= slices; z++) {
            for (int c = 1; c <= channels; c++) {
                ImageProcessor src = source.getStack().getProcessor(source.getStackIndex(c, z, t))
                src.setRoi(x0, y0, w, h)
                ImageProcessor out = src.crop()
                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < w; x++) {
                        int gx = x0 + x, gy = y0 + y, gz = z - 1
                        boolean inNucleus = (int) nucleusLabels.getPixel(gx, gy, gz) == nucleusLabel
                        boolean inSpot = (int) spotLabels.getPixel(gx, gy, gz) == spotLabel
                        boolean keep = c == nucleusChannel ? inNucleus :
                                (c == spotChannel ? inSpot : (inNucleus || inSpot))
                        if (!keep) out.setf(x, y, 0.0f)
                    }
                }
                result.addSlice(out)
                src.resetRoi()
            }
        }
    }
    ImagePlus output = new ImagePlus("masked target crop", result)
    output.setDimensions(channels, slices, frames)
    output.setOpenAsHyperStack(true)
    output.setCalibration(source.getCalibration().copy())
    return output
}

// 某通道的最大值投影（float 数组）。
float[] maxProjection(ImagePlus source, int channel) {
    int w = source.getWidth(); int h = source.getHeight()
    float[] out = new float[w * h]
    Arrays.fill(out, -Float.MAX_VALUE)
    for (int z = 1; z <= source.getNSlices(); z++) {
        ImageProcessor ip = source.getStack().getProcessor(source.getStackIndex(channel, z, 1))
        for (int i = 0; i < w * h; i++) {
            float v = ip.getf(i % w, (int) (i / w))
            if (v > out[i]) out[i] = v
        }
    }
    return out
}

// 某通道在 XY 窗口内、全 z 的最大值投影，直接从原栈读，不复制裁切。
float[] windowMaxProjection(ImagePlus source, int channel, int x0, int y0, int w, int h) {
    float[] out = new float[w * h]
    Arrays.fill(out, -Float.MAX_VALUE)
    for (int z = 1; z <= source.getNSlices(); z++) {
        ImageProcessor ip = source.getStack().getProcessor(source.getStackIndex(channel, z, 1))
        for (int y = 0; y < h; y++) {
            int base = y * w
            for (int x = 0; x < w; x++) {
                float v = ip.getf(x0 + x, y0 + y)
                if (v > out[base + x]) out[base + x] = v
            }
        }
    }
    return out
}

double[] percentiles(float[] values, double low, double high) {
    float[] copy = values.clone()
    Arrays.sort(copy)
    int n = copy.length
    double lo = copy[(int) Math.max(0, Math.min(n - 1, Math.round(low * (n - 1))))]
    double hi = copy[(int) Math.max(0, Math.min(n - 1, Math.round(high * (n - 1))))]
    if (hi <= lo) hi = lo + 1
    return [lo, hi] as double[]
}

int toByte(double value) {
    return (int) Math.max(0, Math.min(255, Math.round(value)))
}

// 把 RGB 的蓝/红两个通道画进画布的一块区域。与 1.1/1.2 显示参数相同。
void paintTile(ColorProcessor canvas, float[] blue, float[] red, int px, int py, int w, int h,
               double clipLow, double clipHigh) {
    double[] blueRange = percentiles(blue, clipLow, clipHigh)
    double[] redRange = percentiles(red, clipLow, clipHigh)
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int i = y * w + x
            int b = toByte(255.0 * (blue[i] - blueRange[0]) / (blueRange[1] - blueRange[0]))
            int r = toByte(255.0 * (red[i] - redRange[0]) / (redRange[1] - redRange[0]))
            canvas.set(px + x, py + y, (r << 16) | (0 << 8) | b)
        }
    }
}

// 原始整视野决策图。只画 1 px 核边界，不填充、无细胞内文字，避免遮盖 DAPI 或中心体。
// 单个细胞的编号、完整保留/剔除原因仍在 archive 的 *_nuclei_all.csv 中查询。
void saveDecisionMap(ImagePlus raw, ImageInt labels, List orderedAll, Map newIdOf,
                     File target, String stem, int pad, double clipLow, double clipHigh,
                     int nucleusChannel, int spotChannel) {
    int w = raw.getWidth(), h = raw.getHeight(), slices = raw.getNSlices()
    int titleH = 54
    ColorProcessor canvas = new ColorProcessor(w + 2 * pad, h + 2 * pad + titleH)
    canvas.setColor(new Color(18, 18, 18))
    canvas.fill()
    canvas.setAntialiasedText(true)
    paintTile(canvas, maxProjection(raw, nucleusChannel), maxProjection(raw, spotChannel),
              pad, titleH + pad, w, h, clipLow, clipHigh)

    int keptCount = 0, droppedCount = 0
    orderedAll.each { Object3D nucleus ->
        int label = nucleus.getValue()
        boolean selected = newIdOf.containsKey(label)
        if (selected) keptCount++; else droppedCount++
        Color color = selected ? new Color(70, 255, 140) : new Color(255, 220, 50)
        boolean[] outline = borderOf(insideProjection(labels, label, 0, 0, w, h, slices), w, h)
        canvas.setColor(color)
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (outline[y * w + x]) canvas.set(pad + x, titleH + pad + y, color.getRGB())
            }
        }
    }

    canvas.setFont(new Font("SansSerif", Font.BOLD, 14))
    canvas.setColor(Color.WHITE)
    drawFit(canvas, "${stem}.nd2 原图选择 / 剔除决策图 (split v2.0)", pad, 18, w)
    canvas.setFont(new Font("SansSerif", Font.PLAIN, 12))
    canvas.setColor(new Color(70, 255, 140))
    drawFit(canvas, "绿色细边 = kept ${keptCount}", pad, 36, (int) (w * 0.42))
    canvas.setColor(new Color(255, 220, 50))
    drawFit(canvas, "黄色细边 = dropped ${droppedCount}", pad + (int) (w * 0.42), 36, (int) (w * 0.34))
    canvas.setColor(new Color(180, 180, 180))
    drawFit(canvas, "原图：DAPI 蓝，γ-tubulin 红；编号/原因见 CSV", pad, 50, w)
    new FileSaver(new ImagePlus("decision map", canvas)).saveAsPng(target.absolutePath)
}

// ---------- 文字排版：一律实测宽度，超宽就截断 ----------
// ImageJ 的 ColorProcessor.getStringWidth() 走的就是当前 Font 的 FontMetrics，
// 所以这里量出来的宽度和实际画出来的一致。1.3 是凭估计排的，结果两格的字互相压住。
String fitString(ColorProcessor cp, String text, int maxWidth) {
    if (text == null || text.isEmpty()) return ""
    if (cp.getStringWidth(text) <= maxWidth) return text
    String t = text
    while (t.length() > 1 && cp.getStringWidth(t + "...") > maxWidth) {
        t = t.substring(0, t.length() - 1)
    }
    TRUNCATED_STRINGS++
    return t + "..."
}

// 画一行，超出 maxWidth 就截断；返回实际画出来的宽度，便于在同一行后面接着画别的。
int drawFit(ColorProcessor cp, String text, int x, int y, int maxWidth) {
    String t = fitString(cp, text, maxWidth)
    cp.drawString(t, x, y)
    return cp.getStringWidth(t)
}

// 网格总览（保留细胞）：每格一个细胞，DAPI 蓝 + 中心体通道红的最大值投影。
// 图像区域一个像素都不覆盖：中心体本身就是要目检的红点，早先版本在上面画过十字
// 标记会正好盖住信号，用户 2026-09-20 要求去掉。
// 1.4：每一行文字都走 drawFit 实测截断，格子最小宽度 200。
void saveOverview(List records, File target, int columns, int pad, int header, int minTileW,
                  double clipLow, double clipHigh) {
    if (records.isEmpty()) return
    int cols = Math.min(columns, records.size())
    int rows = (int) Math.ceil(records.size() / (double) cols)
    int tileW = Math.max(minTileW, records.collect { it.w }.max() + 2 * pad)
    int tileH = records.collect { it.h }.max() + 2 * pad + header
    int textW = tileW - 2 * pad
    ColorProcessor canvas = new ColorProcessor(cols * tileW, rows * tileH)
    canvas.setColor(new Color(18, 18, 18))
    canvas.fill()
    canvas.setAntialiasedText(true)

    records.eachWithIndex { record, int index ->
        int col = index % cols
        int row = (int) (index / cols)
        int ox = col * tileW
        int oy = row * tileH
        int w = record.w; int h = record.h
        int px = ox + pad + (int) ((tileW - 2 * pad - w) / 2)
        int py = oy + header + pad
        paintTile(canvas, record.blue, record.red, px, py, w, h, clipLow, clipHigh)

        boolean contaminated = !record.flags.isEmpty()
        Color frame = contaminated ? new Color(255, 80, 255) :
                (record.split ? new Color(80, 200, 255) :
                (record.fused ? new Color(255, 150, 40) :
                (record.edge ? new Color(255, 230, 60) : new Color(90, 90, 90))))
        canvas.setColor(frame)
        canvas.drawRect(px - 1, py - 1, w + 2, h + 2)

        // 第 1 行：文件名 + 旧编号。两段分开量宽，第二段用第一段剩下的宽度。
        canvas.setFont(new Font("SansSerif", Font.BOLD, 13))
        canvas.setColor(Color.WHITE)
        int used = drawFit(canvas, record.name.replaceFirst(/\.tif$/, ""), ox + pad, oy + pad + 13, textW)
        canvas.setFont(new Font("SansSerif", Font.PLAIN, 11))
        canvas.setColor(new Color(150, 150, 150))
        drawFit(canvas, " (v1.0 #${record.oldId})", ox + pad + used + 4, oy + pad + 13,
                textW - used - 4)
        // 第 2 行：尺寸 + 标记
        def tags = []
        if (record.split) tags << "SPLIT"
        if (record.manual) tags << "MANUAL"
        if (record.edge) tags << "EDGE"
        if (record.fused) tags << "FUSED?"
        canvas.setColor(new Color(170, 170, 170))
        drawFit(canvas, "${w}x${h} px" + (tags.isEmpty() ? "" : "  " + tags.join(" ")),
                ox + pad, oy + pad + 27, textW)
        // 第 3 行：距离 + 亮度
        canvas.setColor(new Color(120, 255, 160))
        drawFit(canvas, "MTOC ${String.format('%.2f', record.distance)} um  maxI ${record.maxI}",
                ox + pad, oy + pad + 40, textW)
        // 第 4 行：污染标记
        if (contaminated) {
            canvas.setColor(frame)
            canvas.setFont(new Font("SansSerif", Font.BOLD, 11))
            drawFit(canvas, record.flags.join(" "), ox + pad, oy + pad + 53, textW)
        }
    }
    new FileSaver(new ImagePlus("overview", canvas)).saveAsPng(target.absolutePath)
}

// 被剔除细胞缩略图：固定窗口，画被剔除核自身的轮廓，外加最多一个候选圆环。
// 1.4 相对 1.3 的三处改动（用户 2026-09-20 目检后要求）：
//   - 标注每项一行，全部实测宽度截断，不再两项挤一行、串到隔壁格；
//   - 只给“该核自己分到的最亮候选”画圆环并标一次归属；该核一个候选都没分到时，
//     改画窗口内最亮的那个并标出它归属给谁。其余候选只画 2 px 小点，不带文字。
//     1.3 在 #19 那格画了 22 个带标签的圆环，把细胞整个盖住了；
//   - 表头说明分两行。
void saveDroppedPreview(List tiles, File target, String stem, int columns, int pad,
                        int markerRadius, int dotRadius, int minTileW,
                        double clipLow, double clipHigh, int gate) {
    if (tiles.isEmpty()) return
    int header = 84
    int cols = Math.min(columns, tiles.size())
    int rows = (int) Math.ceil(tiles.size() / (double) cols)
    int tileW = Math.max(minTileW, tiles.collect { it.w }.max() + 2 * pad)
    int tileH = tiles.collect { it.h }.max() + 2 * pad + header
    int textW = tileW - 2 * pad
    int titleBar = 58
    ColorProcessor canvas = new ColorProcessor(cols * tileW, rows * tileH + titleBar)
    canvas.setColor(new Color(18, 18, 18))
    canvas.fill()
    canvas.setAntialiasedText(true)
    int fullW = cols * tileW - 2 * pad
    canvas.setFont(new Font("SansSerif", Font.BOLD, 14))
    canvas.setColor(new Color(255, 120, 120))
    drawFit(canvas, "${stem}.nd2  split v1.5 剔除细胞 - ${tiles.size()} 个核",
            pad, 18, fullW)
    canvas.setFont(new Font("SansSerif", Font.PLAIN, 12))
    canvas.setColor(new Color(180, 180, 180))
    drawFit(canvas, "青色轮廓 = 被剔除核；圆环 = 最亮候选（双环表示 maxI >= ${gate}）；#n = 候选归属核。", pad, 36, fullW)
    drawFit(canvas, "小点 = 窗口内其他候选（不标号）。每格下方为中文剔除原因；窗口居中于该核。", pad, 52, fullW)

    tiles.eachWithIndex { t, int index ->
        int col = index % cols
        int row = (int) (index / cols)
        int ox = col * tileW
        int oy = row * tileH + titleBar
        int w = t.w, h = t.h
        int px = ox + pad + (int) ((tileW - 2 * pad - w) / 2)
        int py = oy + header + pad
        paintTile(canvas, t.blue, t.red, px, py, w, h, clipLow, clipHigh)

        // 被剔除核自身的轮廓：1 px 青线。
        boolean[] outline = t.outline
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (outline[y * w + x]) canvas.set(px + x, py + y, (0 << 16) | (230 << 8) | 230)
            }
        }
        // 其余候选：2 px 小点，暗黄色，不带文字。
        canvas.setColor(new Color(150, 130, 60))
        t.marks.each { m ->
            if (m.ringed) return
            canvas.fillOval(px + m.x - dotRadius, py + m.y - dotRadius, 2 * dotRadius, 2 * dotRadius)
        }
        // 唯一的圆环：该核自己的最亮候选（没有就是窗口内最亮的那个）。
        def ring = t.marks.find { it.ringed }
        if (ring != null) {
            Color c = ring.mine ? (ring.bright ? new Color(120, 255, 160) : new Color(255, 170, 60))
                                : new Color(255, 120, 255)
            canvas.setColor(c)
            canvas.drawOval(px + ring.x - markerRadius, py + ring.y - markerRadius,
                            2 * markerRadius, 2 * markerRadius)
            if (ring.bright) {
                canvas.drawOval(px + ring.x - markerRadius - 2, py + ring.y - markerRadius - 2,
                                2 * markerRadius + 4, 2 * markerRadius + 4)
            }
            // 标签默认画在圆环右侧；右边装不下就翻到左侧，不要截成 "#..." 把归属编号丢掉。
            canvas.setFont(new Font("SansSerif", Font.BOLD, 11))
            String lab = "#" + ring.owner
            int labW = canvas.getStringWidth(lab)
            int lx = px + ring.x + markerRadius + 3
            if (lx + labW > px + w) lx = px + ring.x - markerRadius - 3 - labW
            if (lx < px) lx = px
            drawFit(canvas, lab, lx, py + ring.y + 4, px + w - lx)
        }

        canvas.setColor(new Color(255, 90, 90))
        canvas.drawRect(px - 1, py - 1, w + 2, h + 2)

        // 标注：四行，每行一项，全部实测截断。
        canvas.setFont(new Font("SansSerif", Font.BOLD, 13))
        canvas.setColor(Color.WHITE)
        int used = drawFit(canvas, "剔除 v1.0 #${t.oldId}", ox + pad, oy + pad + 13, textW)
        if (t.edge) {
            canvas.setFont(new Font("SansSerif", Font.BOLD, 11))
            canvas.setColor(new Color(255, 230, 60))
            drawFit(canvas, " EDGE", ox + pad + used + 4, oy + pad + 13, textW - used - 4)
        }
        canvas.setFont(new Font("SansSerif", Font.PLAIN, 11))
        canvas.setColor(new Color(170, 170, 170))
        drawFit(canvas, "核体积 ${String.format('%.0f', t.volUm3)} um3", ox + pad, oy + pad + 29, textW)
        canvas.setColor(new Color(120, 255, 160))
        drawFit(canvas, t.maxI < 0 ? "本核无中心体候选" :
                "本核最亮候选 maxI ${t.maxI}", ox + pad, oy + pad + 43, textW)
        canvas.setColor(new Color(150, 200, 255))
        drawFit(canvas, Double.isNaN(t.near) ? "最近候选：无" :
                "最近候选 ${String.format('%.2f', t.near)} um", ox + pad, oy + pad + 57, textW)
        canvas.setFont(new Font("SansSerif", Font.BOLD, 11))
        canvas.setColor(new Color(255, 120, 120))
        drawFit(canvas, shortReason(t.reason), ox + pad, oy + pad + 71, textW)
    }
    new FileSaver(new ImagePlus("dropped", canvas)).saveAsPng(target.absolutePath)
}

// 递归拆分对照图：每个候选对象一行；一行里第 0 格是拆分前，第 k 格是第 k 轮之后。
// 块的颜色按第一次出现的顺序固定，所以同一块在相邻两格里是同色，被切开时才换成两个新色。
void saveSplitQc(List tiles, File target, String stem, int pad, double clipLow, double clipHigh,
                 int medianVox, double ratio, float tolUm) {
    if (tiles.isEmpty()) return
    def partColors = [new Color(80, 255, 140), new Color(255, 170, 40), new Color(255, 90, 220),
                      new Color(90, 190, 255), new Color(255, 255, 90), new Color(160, 130, 255),
                      new Color(120, 255, 255), new Color(255, 130, 110)]
    int maxPanels = Math.max(2, tiles.collect { it.panels.size() }.max())
    int header = 84
    int titleBar = 58
    int tileW = Math.max(300, tiles.collect { it.w }.max() + 2 * pad)
    int tileH = tiles.collect { it.h }.max() + 2 * pad + header
    int textW = tileW - 2 * pad
    ColorProcessor canvas = new ColorProcessor(maxPanels * tileW, tiles.size() * tileH + titleBar)
    canvas.setColor(new Color(18, 18, 18))
    canvas.fill()
    canvas.setAntialiasedText(true)
    int fullW = maxPanels * tileW - 2 * pad
    canvas.setFont(new Font("SansSerif", Font.BOLD, 14))
    canvas.setColor(new Color(120, 220, 255))
    drawFit(canvas, "${stem}.nd2  recursive 3D watershed split QC (v2.0) - ${tiles.size()} object(s) above " +
            "${ratio} x median (${medianVox} vox)", pad, 18, fullW)
    canvas.setFont(new Font("SansSerif", Font.PLAIN, 12))
    canvas.setColor(new Color(180, 180, 180))
    drawFit(canvas, "column 0 = before any split; column k = after round k.  one colour per piece, " +
            "kept across columns, so a piece that gets cut shows up as two new colours.", pad, 36, fullW)
    drawFit(canvas, "outline + 35% tint per piece.  EDT + MaximaFinder(radXY 2 px, radZ 5 px, tol ${tolUm} um) " +
            "+ Watershed3D, re-applied to any piece still above the cutoff.", pad, 52, fullW)

    tiles.eachWithIndex { t, int index ->
        def rec = t.rec
        int oy = index * tileH + titleBar
        int w = t.w, h = t.h

        for (int col = 0; col < maxPanels; col++) {
            int ox = col * tileW
            int px = ox + pad + (int) ((tileW - 2 * pad - w) / 2)
            int py = oy + header + pad
            // 超出该对象实际轮数的列留空，只写一句说明。
            if (col >= t.panels.size()) {
                canvas.setFont(new Font("SansSerif", Font.PLAIN, 11))
                canvas.setColor(new Color(110, 110, 110))
                drawFit(canvas, "(no round ${col} for this object)", ox + pad, oy + pad + 13, textW)
                continue
            }
            def panel = t.panels[col]
            paintTile(canvas, t.blue, t.red, px, py, w, h, clipLow, clipHigh)
            panel.each { p ->
                Color c = partColors[p.colorIx % partColors.size()]
                boolean[] inside = p.inside
                boolean[] border = borderOf(inside, w, h)
                for (int i = 0; i < w * h; i++) {
                    int cx = px + (i % w), cy = py + (int) (i / w)
                    if (border[i]) {
                        canvas.set(cx, cy, (c.getRed() << 16) | (c.getGreen() << 8) | c.getBlue())
                    } else if (inside[i]) {
                        int v = canvas.get(cx, cy)
                        int r = (int) (0.65 * ((v >> 16) & 255) + 0.35 * c.getRed())
                        int g = (int) (0.65 * ((v >> 8) & 255) + 0.35 * c.getGreen())
                        int b = (int) (0.65 * (v & 255) + 0.35 * c.getBlue())
                        canvas.set(cx, cy, (r << 16) | (g << 8) | b)
                    }
                }
            }
            canvas.setColor(new Color(90, 90, 90))
            canvas.drawRect(px - 1, py - 1, w + 2, h + 2)

            canvas.setFont(new Font("SansSerif", Font.BOLD, 13))
            canvas.setColor(Color.WHITE)
            drawFit(canvas, col == 0 ? "BEFORE  v1.0 #${rec.oldId}" : "AFTER ROUND ${col}",
                    ox + pad, oy + pad + 13, textW)
            canvas.setFont(new Font("SansSerif", Font.PLAIN, 11))
            canvas.setColor(new Color(170, 170, 170))
            drawFit(canvas, "${panel.size()} piece(s), cutoff ${(int) (ratio * medianVox)} vox",
                    ox + pad, oy + pad + 28, textW)
            // 每块一行：编号、体素、占母对象比例、是否仍然超标。
            int line = 0
            panel.sort(false) { -it.voxels }.each { p ->
                if (line >= 4) return
                canvas.setColor(partColors[p.colorIx % partColors.size()])
                drawFit(canvas, "#${p.id}: ${p.voxels} vox " +
                        "(${String.format('%.0f', 100.0 * p.voxels / rec.parentVox)}%, " +
                        "${String.format('%.2f', p.voxels / (double) medianVox)}x med)" +
                        (p.voxels > ratio * medianVox ? "  STILL OVER" : ""),
                        ox + pad, oy + pad + 42 + 13 * line, textW)
                line++
            }
            if (panel.size() > 4) {
                canvas.setColor(new Color(150, 150, 150))
                drawFit(canvas, "... ${panel.size() - 4} more", ox + pad, oy + pad + 42 + 13 * 4, textW)
            }
        }
    }
    new FileSaver(new ImagePlus("split qc", canvas)).saveAsPng(target.absolutePath)
}

// 剔除理由缩成一行，长句会压到隔壁格上。
String shortReason(String reason) {
    if (reason == null || reason.trim().isEmpty()) return ""
    def parts = []
    if (reason.contains("未分配到本核")) parts << "未分到本核中心体"
    if (reason.contains("强度")) parts << "中心体亮度不足"
    if (reason.contains("距离")) parts << "中心体距离过远"
    if (reason.contains("疑似粘连核")) parts << "核仍疑似粘连"
    if (reason.contains("邻近中心体")) parts << "裁切含邻中心体"
    if (reason.contains("邻近细胞核")) parts << "裁切含邻核"
    return parts.isEmpty() ? reason : parts.join("；")
}
