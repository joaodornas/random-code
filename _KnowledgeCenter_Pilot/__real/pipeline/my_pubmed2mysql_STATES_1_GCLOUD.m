function my_pubmed2mysql_STATES_1_GCLOUD

query_id = 1;

start_year = 1900;
end_year = 2016;

I_Am_Testing = 0;
Environment = 'GCLOUD';

start_year = 1900;
end_year = 2016;
query_label = 'awake';
my_pubmed2mysql_v2(query_label,query_id,start_year,end_year,I_Am_Testing,Environment);

start_year = 1900;
end_year = 2016;
query_label = 'resting_state';
my_pubmed2mysql_v2(query_label,query_id,start_year,end_year,I_Am_Testing,Environment);

start_year = 1900;
end_year = 2016;
query_label = 'sleep';
my_pubmed2mysql_v2(query_label,query_id,start_year,end_year,I_Am_Testing,Environment);


end

