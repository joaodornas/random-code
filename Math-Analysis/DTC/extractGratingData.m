function extractGratingData

nConditions = 16;

DTC_size = 200;

lambda(1) = 1;
lambda(2) = 2;
lambda(3) = 4;
lambda(4) = 8;
lambda(5) = 10;

for l=1:5
    
    gratingAnyScale(lambda(l),DTC_size,nConditions);
    
end

function gratingAnyScale(scale,DTC_size,nConditions)

        for i=-2:2

            for j=-1:2
            
                if (i ~= 2) && (j ~= -1)
                    
                    DTC_folderpath = strcat('/DTC-200_66cm_Sf',int2str(i),',0_Tf',int2str(j),',0_C0,95_L1,0/');

                    fullpath = strcat('/Volumes/Data/DATA/DTC/protocols',DTC_folderpath);

                    for n=1:nConditions

                        if n < 10

                            number = strcat('0',int2str(n));

                        else

                            number = int2str(n);

                        end

                        filepath = strcat(fullpath,'Condition',number,'-0.bmp');

                        img = imread(filepath);

                        grating = img(1:DTC_size,1:DTC_size);

                        grating_scaled = imresize(grating,[DTC_size/scale DTC_size/scale],'bicubic');

                        gratingFourier = fft2(grating_scaled);

                        scale_path = strcat('scale-',int2str(scale),'/');

                        fourier_path = strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path);

                        mkdir(strcat(fourier_path,'Sf',int2str(i),'_Tf',int2str(j)));

                        grating_path = strcat('Sf',int2str(i),'_Tf',int2str(j));

                        imwrite(grating,strcat('/Volumes/Data/DATA/DTC/Gratings/',scale_path,'Grating-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n),'.jpg'),'jpeg');

                        imwrite(grating_scaled,strcat('/Volumes/Data/DATA/DTC/Gratings/',scale_path,'Grating_Scaled-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n),'.jpg'),'jpeg');

                        imwrite(gratingFourier,strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Fourier-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n),'.jpg'),'jpeg');

                        save(strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Fourier-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n),'.mat'),'gratingFourier');

                        A = figure;

                        amplitudes = abs(gratingFourier);

                        mesh(amplitudes);

                        print(A,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Amplitudes-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        close all;

                        G = figure;

                        angulos = angle(gratingFourier);

                        mesh(angulos);

                        print(G,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Angulos-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        close all;

                        angulos_negativos = getNegativos(angulos);

                        angulos_positivos = getPositivos(angulos);

                        angulos_nulos = getNulos(angulos);

                        name = strcat('Grating-','SF-',int2str(i),'-TF-',int2str(j));

                        fourier_info_struct = struct('name',name,'gratingFourier',gratingFourier,'amplitudes',amplitudes,'angulos',angulos,'angulos_negativos',angulos_negativos,'angulos_positivos',angulos_positivos,'angulos_nulos',angulos_nulos);

                        save(strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Fourier-Info-Struct-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n),'.mat'),'fourier_info_struct');

                        G_negativos = figure;

                        mesh(angulos_negativos);

                        print(G_negativos,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Mesh-Angulos-Negativos-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        G_positivos = figure;

                        mesh(angulos_positivos);

                        print(G_positivos,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Mesh-Angulos-Positivos-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        G_nulos = figure;

                        mesh(angulos_nulos);

                        print(G_nulos,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Mesh-Angulos-Nulos-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        close all;

                        f = figure;

                        image(grating);
                        colormap(gray); 
                        hold on;

                        dimension_X = size(angulos_negativos,1);
                        dimension_Y = size(angulos_negativos,2);

                        for s=1:dimension_X

                            for t=1:dimension_Y

                                handle = text(s*scale,t*scale,num2str(angulos_negativos(s,t)));
                                set(handle,'Color','r');

                            end

                        end

                        print(f,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Sobre-Angulos-Negativos-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        g = figure;

                        image(grating);
                        colormap(gray);
                        hold on;

                        dimension_X = size(angulos_positivos,1);
                        dimension_Y = size(angulos_positivos,2);

                        for s=1:dimension_X

                            for t=1:dimension_Y

                                handle = text(s*scale,t*scale,num2str(angulos_positivos(s,t)));
                                set(handle,'Color','b');

                            end

                        end
                        print(g,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Sobre-Angulos-Positivos-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        h = figure;

                        image(grating);
                        colormap(gray);
                        hold on;

                        dimension_X = size(angulos_nulos,1);
                        dimension_Y = size(angulos_nulos,2);

                        for s=1:dimension_X

                            for t=1:dimension_Y

                                handle = text(s*scale,t*scale,num2str(angulos_nulos(s,t)));
                                set(handle,'Color','g');

                            end

                        end
                        print(h,'-depsc',strcat('/Volumes/Data/DATA/DTC/Gratings-Fourier/',scale_path,grating_path,'/','Grating-Sobre-Angulos-Nulos-SF-',int2str(i),'-TF-',int2str(j),'-C0,95_L1,0-condition-',int2str(n)));

                        close all;
                        
                    end
                
                end
                
            end

        end

    end

end



function newMatrix = getNegativos(matrix)

    dimension_1 = size(matrix,1);
    dimension_2 = size(matrix,2);
    
    for x=1:dimension_1
        
        for y=1:dimension_2
            
            if matrix(x,y) >= 0
                
                matrix(x,y) = 0;
                
            end
            
        end
        
    end

    
    newMatrix = matrix;
    
end

function newMatrix = getPositivos(matrix)

    dimension_1 = size(matrix,1);
    dimension_2 = size(matrix,2);
    
    for x=1:dimension_1
        
        for y=1:dimension_2
            
            if matrix(x,y) <= 0
                
                matrix(x,y) = 0;
                
            end
            
        end
        
    end

    
    newMatrix = matrix;
    
end

function newMatrix = getNulos(matrix)

    dimension_1 = size(matrix,1);
    dimension_2 = size(matrix,2);
    
    for x=1:dimension_1
        
        for y=1:dimension_2
            
            if matrix(x,y) ~= 0
                
                matrix(x,y) = 0;
                
            else
                
                matrix(x,y) = 100;
                
            end
            
        end
        
    end

    
    newMatrix = matrix;
    
end



