function sc=speckle_contrast(v)
x=double(v);
sc=std(x)./mean(x);
end
