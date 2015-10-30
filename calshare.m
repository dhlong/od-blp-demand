function s = calshare(delta, emu, iT)

emu = bsxfun(@times, emu, exp(delta));
diT = diff(iT);

if all(diT<=1) && all(diT>=0)
    % if iT is in running order    
    cdindex = [find(diff(iT)>0); length(iT)];
    s0 = cumsum(emu);
    s0 = diff([zeros(size(s0(1,:)));s0(cdindex,:)]) + 1;
else
    % if iT is not in running order
    [ii, jj] = meshgrid(1:size(emu,2), iT);
    ii = ii(:); jj = jj(:);
    s0 = 1 + accumarray([jj ii], emu(:));    
end

s = emu./s0(iT,:);

end