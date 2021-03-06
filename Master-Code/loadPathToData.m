function [pathToData,dataFolders] = loadPathToData(protocol)

   
    switch protocol
        
        case 'VideosFB'
            
            pathToData = 'C:\Users\Jo�o V. Dornas\Hard Quale .com\DornasLab.org - Documents\_CODE\TP-RANDOM\Master-Code\DATA\FB\Registro\';
     
            dataFolders = {'SUA - CELL 1 - 12-08-20_\nsp003a01\_nsp003a01_1b-v3.mat'...
                'SUA - CELL 2 - 12-08-27_\nsp005a01\_nsp005a01_2b-v3.mat'...
                'SUA - CELL 3 - 12-08-29_\nsp006a01\_nsp006a01_1b-v3.mat'...
                'SUA - CELL 4 - 12-08-29_\nsp006b01\_nsp006b01_1b-v1.mat'...
                'SUA - CELL 5 - 12-08-30_\nsp007a01\_nsp007a01_1b-v4.mat'...
                'SUA - CELL 5 - 12-08-30_\nsp007a02\_nsp007a02_1b-v1.mat'...
                'SUA - CELL 6 - 12-08-30_\nsp007a01\_nsp007a01_2b-v4.mat'...
                'SUA - CELL 6 - 12-08-30_\nsp007a02\_nsp007a02_2b-v1.mat'...
                'SUA - CELL 6 - 12-08-30_\nsp007a03\_nsp007a03_2a-v3.mat'...
                'SUA - CELL 7 - 12-08-31_\nsp008a01\_nsp008a01_1b-v1.mat'...
                'SUA - CELL 7 - 12-08-31_\nsp008a02\_nsp008a02_1a-v3.mat'...
                'SUA - CELL 7 - 12-08-31_\nsp008a03\_nsp008a03_1b-v4.mat'...
                'SUA - CELL 8 - 12-08-31_\nsp008a01\_nsp008a01_2b-v1.mat'...
                'SUA - CELL 8 - 12-08-31_\nsp008a02\_nsp008a02_2b-v3.mat'...
                'SUA - CELL 8 - 12-08-31_\nsp008a03\_nsp008a03_2a-v4.mat'...
                'SUA - CELL 9 - 12-09-03_\nsp009a01\_nsp009a01_2b-v2.mat'...
                'SUA - CELL 10 - 12-09-03_\nsp009b01\_nsp009b01_1a-v4.mat'...
                'SUA - CELL 11 - 12-09-04_\nsp010a01\_nsp010a1_1b-v2.mat'...
                'SUA - CELL 11 - 12-09-04_\nsp010a02\_nsp010a02_1b-v4.mat'...
                'SUA - CELL 12 - 12-09-05_\nsp011a01\_nsp011a01_1a-v3.mat'...
                'SUA - CELL 12 - 12-09-05_\nsp011a02\_nsp011a02_1b-v4.mat'...
                'SUA - CELL 13 - 12-09-06_\nsp012a01\_nsp012a01_1b-v3.mat'...
                'SUA - CELL 13 - 12-09-06_\nsp012a02\_nsp012a02_1a-v4.mat'...
                'SUA - CELL 14 - 12-10-16_\nsp013a01\_nps013a01_1b-v1.mat'...
                'SUA - CELL 15 - 12-12-17_\nsp033a01\_nsp033a01_1b.mat'...
                'SUA - CELL 15 - 12-12-17_\nsp033a09\_nsp033a09_1b-v311.mat'...
                'SUA - CELL 16 - 12-12-17_\nsp033a01\_nsp033a01_1c.mat'...
                'SUA - CELL 16 - 12-12-17_\nsp033a09\_nsp033a09_1c-v321.mat'...
                'SUA - CELL 17 - 12-12-17_\nsp033a01\_nsp033a01_3a.mat'...
                'SUA - CELL 17 - 12-12-17_\nsp033a09\_nsp033a09_3b-v331.mat'};
            
        case 'VideosCC'
            
             pathToData = 'C:\Users\Jo�o V. Dornas\Hard Quale .com\DornasLab.org - Documents\_CODE\TP-RANDOM\Master-Code\DATA\CC\Registro\';

             dataFolders = {'SUA - CELL 1 - 12-10-22_\nsp015a01\_nsp015a01_2a-v1.mat'...
                'SUA - CELL 2 - 12-10-22_\nsp015a02\_nsp015a02_2b-v3.mat'...
                'SUA - CELL 3 - 12-10-23_\nsp016a01\_nsp016a01_2b-v2.mat'...
                'SUA - CELL 4 - 12-10-24_\nsp016a02\_nsp016a02_2b-v4.mat'...
                'SUA - CELL 5 - 12-10-24_\nsp017b01\_nsp017b01_2b-v2.mat'...
                'SUA - CELL 6 - 12-10-26_\nsp017a01\_nsp017a01_1a-v3.mat'...
                'SUA - CELL 6 - 12-10-26_\nsp018b01\_nsp018b01_1b-v3.mat'...
                'SUA - CELL 6 - 12-10-26_\nsp018b02\_nsp018b02_1b.mat'...
                'SUA - CELL 7 - 12-11-01_\nsp020a01\_nsp020a01_1b-v5.mat'...
                'SUA - CELL 7 - 12-11-01_\nsp020a02\_nsp020a02_1b-v6.mat'...
                'SUA - CELL 7 - 12-11-01_\nsp020a03\_nsp020a03_1b.mat'...
                'SUA - CELL 8 - 12-11-06_\nsp021a01\_nsp021a01_1a-v5.mat'...
                'SUA - CELL 8 - 12-11-06_\nsp021a02\_nsp021a02_1a-v6.mat'...
                'SUA - CELL 8 - 12-11-06_\nsp021a03\_nsp021a03_1a.mat'...
                'SUA - CELL 9 - 12-11-07_\nsp022a01\_nsp022a01_1b-v5.mat'...
                'SUA - CELL 9 - 12-11-07_\nsp022a02\_nsp022a02_1b-v6.mat'...
                'SUA - CELL 9 - 12-11-07_\nsp022a03\_nsp022a03_1b.mat'...
                'SUA - CELL 10 - 12-11-08_\nsp023a01\_nsp023a01_1a-v1.mat'...
                'SUA - CELL 10 - 12-11-08_\nsp023a02\_nsp023a02_1a-v4'...
                'SUA - CELL 10 - 12-11-08_\nsp023a03\_nsp023a03_1b.mat'...
                'SUA - CELL 13 - 12-12-03_\nsp026a01\_nsp026a01_3b-v8'...
                'SUA - CELL 13 - 12-12-03_\nsp026a02\_nsp026a02_3a-v11'...
                'SUA - CELL 13 - 12-12-03_\nsp026a03\_nsp026a03_3b.mat'...
                'SUA - CELL 14 - 12-12-04_\nsp027a01\_nsp027a01_1b.mat'...
                'SUA - CELL 14 - 12-12-04_\nsp027a02\_nsp027a02_1b-v0.mat'...
                'SUA - CELL 14 - 12-12-04_\nsp027a03\_nsp027a03_1a-v10.mat'...
                'SUA - CELL 14 - 12-12-04_\nsp027a04\_nsp027a04_1b-v9.mat'...
                'SUA - CELL 15 - 12-12-05_\nsp028a01\_nsp028a01_1a.mat'...
                'SUA - CELL 15 - 12-12-05_\nsp028a02\_nsp028a02_1a-v0.mat'...
                'SUA - CELL 15 - 12-12-05_\nsp028a03\_nsp028a03_1a-v6.mat'...
                'SUA - CELL 15 - 12-12-05_\nsp028a04\_nsp028a04_1b-v8.mat'...
                'SUA - CELL 16 - 12-12-07_\nsp029a01\_nsp029a01_1b.mat'...
                'SUA - CELL 16 - 12-12-07_\nsp029a03\_nsp029a03_1a-v8.mat'...
                'SUA - CELL 16 - 12-12-07_\nsp029a04\_nsp029a04_1b-v3.mat'...
                'SUA - CELL 17 - 12-12-10_\nsp030a01\_nsp030a01_2b.mat'...
                'SUA - CELL 17 - 12-12-10_\nsp030a03\_nsp030a03_2a-v5.mat'...
                'SUA - CELL 17 - 12-12-10_\nsp030a04\_nsp030a04_2b-v8.mat'...
                'SUA - CELL 18 - 12-12-12_\nsp031b01\_nsp031b01_1b.mat'...
                'SUA - CELL 18 - 12-12-12_\nsp031b03\_nsp031b03_1b-v3.mat'...
                'SUA - CELL 18 - 12-12-12_\nsp031b04\_nsp031b04_1b-v4.mat'...
                'SUA - CELL 19 - 12-12-13_\nsp032a01\_nsp032a01_1b.mat'...
                'SUA - CELL 19 - 12-12-13_\nsp032a03\_nsp032a03_1a-v8.mat'...
                'SUA - CELL 19 - 12-12-13_\nsp032a04\_nsp032a04_1a-v2.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a01\_nsp033a01_1b.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a01\_nsp033a01_1c.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a01\_nsp033a01_3a.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a04\_nsp033a04_1b-v821.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a04\_nsp033a04_1c-v811.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a04\_nsp033a04_3a-v831.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a05\_nsp033a05_1b-v511.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a05\_nsp033a05_1c-v521.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a05\_nsp033a05_3a-v531.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a06\_nsp033a06_1b-v812.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a06\_nsp033a06_1c-v822.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a06\_nsp033a06_3b-v832.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a07\_nsp033a07_1b-v512.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a07\_nsp033a07_1c-v522.mat'...
                'SUA - CELL 20_21_22 - 12-12-17_\nsp033a07\_nsp033a07_3b-v532.mat'...
                'SUA - CELL 23_24 - 12-12-18_\nsp034a01\_nsp034a01_2b.mat'...
                'SUA - CELL 23_24 - 12-12-18_\nsp034a04\_nsp034a04_2b-v8.mat'...
                'SUA - CELL 23_24 - 12-12-18_\nsp034a05\_nsp034a05_2b-v11.mat'...
                'SUA - CELL 25 - 12-12-19_\nsp035a01\_nsp035a01_2b.mat'...
                'SUA - CELL 25 - 12-12-19_\nsp035a03\_nsp035a03_2b-v8.mat'...
                'SUA - CELL 25 - 12-12-19_\nsp035a04\_nsp035a04_2b-v1.mat'};
            
    end
        
end