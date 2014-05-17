function check_cpindex(numconparam,cpindex);

if any(cpindex<1) | any(cpindex>numconparam) | any(cpindex~=ceil(cpindex))
  error('Not valid cpindex');
end

