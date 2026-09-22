# Methods overview

The workflow is designed for multichannel 3D fluorescence stacks. A nuclear channel is segmented in 3D, punctate objects are detected in a separate channel, and the selected object is assigned to the nearest nuclear object. The primary distance is calibrated in micrometres from the selected object centroid to the nuclear surface. Peak-to-surface and surface-to-surface diagnostics may be retained.

The code does not silently remove observations because a QC flag is present. Any project-specific threshold or inclusion rule must be declared in the versioned template and described in the methods record.
