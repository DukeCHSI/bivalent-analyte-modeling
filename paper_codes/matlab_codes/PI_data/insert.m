function a = insert(a, b, idx)
a = [a(1:length(a) < idx), b, a(1:length(a) >= idx)];
end