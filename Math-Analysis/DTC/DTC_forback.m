function [start_time, end_time, which_one] = DTC_forback(registro)


switch registro
    
    case '_nsp003a01_1b-v3.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp005a01_1b-v3.mat'
        
        which_one = '_blkdtc002a01_2b.mat';
        start_time = 1000;
        end_time = 11000;
        
    case '_nsp005a01_2b-v3.mat'
        
        which_one = '_blkdtc002a01_2b.mat';
        start_time = 1000;
        end_time = 11000;
        
    case '_nsp006a01_1b-v3.mat' %?
        
        which_one = '_blkdtc003a01_1a.mat';
        start_time = 1000;
        end_time = 8000;

    case '_nsp006a02_1b-v4.mat' %?
        
        which_one = '_blkdtc003a01_1a.mat';
        start_time = 1000;
        end_time = 8000;

    case '_nsp006b01_1b-v1.mat'

        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp007a01_1b-v4.mat'
        
        which_one = '_blkdtc004a01_2b.mat';
        start_time = 1000;
        end_time = 5000;        
        
    case '_nsp007a01_2b-v4.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp007a02_1b-v1.mat'
        
        which_one = '_blkdtc004a01_2b.mat';
        start_time = 1000;
        end_time = 8000;        
        
    case '_nsp007a02_2b-v1.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp007a03_2a-v3.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp008a01_1b-v1.mat'
        
        which_one = '_blkdtc005b01_1b.mat';
        start_time = 500;
        end_time = 2500;        
        
    case '_nsp008a02_1a-v3.mat'
        
        which_one = '_blkdtc005b01_1b.mat';
        start_time = 500;
        end_time = 2500;   
        
    case '_nsp008a03_1b-v4.mat'
        
        which_one = '_blkdtc005b01_1b.mat';
        start_time = 500;
        end_time = 2500;   
        
    case '_nsp008a01_2b-v1.mat'
        
        which_one = '_blkdtc005b01_1b.mat';
        start_time = 500;
        end_time = 2500;   
        
    case '_nsp008a02_2b-v3.mat'
        
        which_one = '_blkdtc005b01_1b.mat';
        start_time = 500;
        end_time = 2500;   
        
    case '_nsp008a03_2a-v4.mat'
        
        which_one = '_blkdtc005b01_1b.mat';
        start_time = 500;
        end_time = 2500;   
        
    case '_nsp009a01_2b-v2.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp009b01_1a-v4.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp010a1_1b-v2.mat' %?
        
        which_one = '_blkdtc006a01_1c.mat';
        start_time = 500;
        end_time = 5000; 
        
    case '_nsp010a02_1b-v4.mat' %?
        
        which_one = '_blkdtc006a01_1c.mat';
        start_time = 500;
        end_time = 5000; 
        
    case '_nsp011a01_1a-v3.mat'
        
        which_one = '_blkdtc007a01_1a.mat';
        start_time = 500;
        end_time = 4500;
        
    case '_nsp011a02_1b-v4.mat'
        
        which_one = '_blkdtc007a01_1a.mat';
        start_time = 500;
        end_time = 4500;
       
    case '_nsp012a01_1b-v3.mat'
        
       which_one = '_blkdtc008a01_1a.mat';
       start_time = 1000;
       end_time = 6000;
       
    case '_nsp012a02_1a-v4.mat'
        
       which_one = '_blkdtc008a01_1a.mat';
       start_time = 500;
       end_time = 4500;
       
    case '_nps013a01_1b-v1.mat'
        
       which_one = '_plastic011a01_1b.mat';
       start_time = 500;
       end_time = 2500;
        
    case '_nsp033a09_1b-v311.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp033a09_1c-v321.mat'
        
        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
    case '_nsp033a09_3b-v331.mat'

        which_one = 'none';
        start_time = 0;
        end_time = 0;
        
end




end

