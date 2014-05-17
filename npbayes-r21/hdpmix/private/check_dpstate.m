function check_dpstate(ppindex,dpstate);

ACTIVE  = 2;
FROZEN  = 1;
HELDOUT = 0;
activeset  = find(dpstate==ACTIVE);
frozenset  = find(dpstate==FROZEN);
heldoutset = find(dpstate==HELDOUT);

if any(sort([activeset frozenset heldoutset])~=1:length(ppindex))
  error('dpstate not active, frozen or heldout');
end
pp      = ppindex(frozenset);
pp      = pp(find(pp>0));
ppstate = dpstate(pp);
if any(ppstate==HELDOUT)
  error('Parent of frozen DP is heldout');
end
pp      = ppindex(activeset);
pp      = pp(find(pp>0));
ppstate = dpstate(pp);
if any(ppstate~=ACTIVE)
  error('Parent of active DP is not active');
end

