w = importdata('1.dat');
tmp = interp1(w(:,1), w(:,2), (0:0.001:1), 'spline');
w = [transpose(0:0.001:1),transpose(tmp)];
y = w(:,1);
uruz = w(:,2);

[val,pos] = max(uruz);
pos
for k = 1:size(y)
    if y(k) > 1.56173585128009e-001
        pos = k;
        break
    end
end
pos

cf_i = trapz(y(1:pos),(1-y(1:pos)).*uruz(1:pos));
cf_o = trapz(y(pos:end),(1-y(pos:end)).*uruz(pos:end));
cf_o/cf_i

