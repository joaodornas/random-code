function [do_I_take_snapshot, snapshot_label] = getRenderingSnapshot

global takeRenderingSnapshot
global current_snapshot_label

do_I_take_snapshot = takeRenderingSnapshot;
snapshot_label = current_snapshot_label;

end

