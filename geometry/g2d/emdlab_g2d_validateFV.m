function emdlab_g2d_validateFV(f,v)

f = f(:);
if any(f<0) || max(f)>size(v,1)
    error('Improper facet-vertices data for geometry.');
end

end
