function [PFDs] = set_PFDs(cells)
%Sets PFDs, in HD space, of a set of N=cells cells
increment = 360.0/cells;
PFDs = 1:cells;
PFDs = PFDs*increment;
end

