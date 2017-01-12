function []=fit_plot(fit_name,bound)

fit_image=fitsread(fit_name);
imshow(fit_image,[0 bound]);

end