% Universidade Federal de Minas Gerais -- UFMG
% Instituto de Ciências Exatas -- ICEx
% Departamento de Ciência da Computação -- DCC
%
% Disciplina: Visão Computacional
% Professor:  Mario Fernando Montenegro Campos

function [lambda1, lambda2, v1] = corners(I, N)
%CORNERS Calcula os autovalores e o autovetor principal da vizinhança dos
%pixels de uma imagem.
%  [lambda1, lambda2, v1] = CORNERS(I, N) calcula os autovalores e o autovetor
%  principal da vizinhança de cada pixel da imagem I. N >= 0 é o tamanho da
%  vizinhança, em pixels (N = 0 considera somente o pixel sem os vizinhos).
%  lambda1 e lambda2 são os autovalores da vizinhança analisada para cada
%  pixel, sendo lambda1 >= lambda2.  v1 contém os autovetores em cada pixel,
%  sendo [v1(i,j,1), v1(i,j,2)] o autovetor calculado para as coordenadas i,j.
%
%  Revisado em 2007/04/16 por Vilar Camara Neto
%  Revisado em 2001/04/02 por José Luiz de S. Pio

I = im2double(I);
%I = rgb2gray(I);
[gr, gc] = gradient(I);
winsize = 2 * N + 1;

% calcula as entradas grr, gcc, grc da matriz simetrica G
onerow = ones(1, winsize) / winsize;
grr = conv2(conv2(gr .^ 2, onerow, 'same'), onerow', 'same');
gcc = conv2(conv2(gc .^ 2, onerow, 'same'), onerow', 'same');
grc = conv2(conv2(gr .* gc, onerow, 'same'), onerow', 'same');

[rows, cols] = size(I);

% traco e determinante das matrizes
tr = grr + gcc;
dt = grr .* gcc - grc .^ 2;

% autovalores e o primeiro autovetor
lambda1 = 0.5 * (tr + sqrt(abs(tr .^2 - 4 * dt)));
lambda2 = 0.5 * (tr - sqrt(abs(tr .^2 - 4 * dt)));
v1 = zeros(size(I, 1), size(I, 2), 2);
v1(:, :, 1) = grr - lambda2;
v1(:, :, 2) = grc;
nrm = sqrt(abs(v1(:, :, 1) .^ 2 + v1(:, :, 2) .^ 2));
nrm(find(nrm == 0)) = 1;
v1(:, :, 1) = v1(:, :, 1) ./ nrm;
v1(:, :, 2) = v1(:, :, 2) ./ nrm;
