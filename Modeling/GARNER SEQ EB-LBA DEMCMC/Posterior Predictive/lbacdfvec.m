function Fj = lbacdfvec(t, v, b, a, s, zlim)

Fj = ones(numel(t), 1);
z = (b - a - t .* v)./(t * s);
z(z > zlim) = zlim;
z(z < -zlim) = -zlim;

p = normcdf(z);
lik = normpdf(z);

Fj = Fj + ((b - a - t .* v)./a) .* p + ((t * s)./a) .* lik;

z = (b - t .* v)./(t .* s);
z(z > zlim) = zlim;
z(z < -zlim) = -zlim;

p = normcdf(z);
lik = normpdf(z);

Fj = Fj - ((b - t .* v)./a) .* p - ((t * s)./a) .* lik;

Fj(Fj > .9999) = .9999;
Fj(Fj < .0001) = .0001;