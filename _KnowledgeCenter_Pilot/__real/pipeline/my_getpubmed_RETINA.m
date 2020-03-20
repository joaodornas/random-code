
function my_getpubmed_RETINA

query_id = 1;
start_year = 1900;
end_year = 2016;
minutes_to_pause = 10;

label = 'retina';
   
    my_getpubmed(label,query_id,start_year,end_year,minutes_to_pause);
    
label = 'retinal';

    my_getpubmed(label,query_id,start_year,end_year,minutes_to_pause);

end