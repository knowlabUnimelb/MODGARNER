function fi = lbapdfvec(t, v, b, a, s, zlim)

fi = 0;
z = (b - a - t .* v)./(t .* s);
z(z > zlim) = zlim;
z(z < -zlim) = -zlim;

p = normcdf(z);
lik = normpdf(z);

fi = fi - v .* p + s .* lik;

z = (b - t .* v)./(t .* s);

z(z > zlim) = zlim;
z(z < -zlim) = -zlim;

p = normcdf(z);
lik = normpdf(z);

fi = fi + v .* p - s .* lik;
fi = (1./a) .* fi;

fi(fi < .0001) = .0001;
