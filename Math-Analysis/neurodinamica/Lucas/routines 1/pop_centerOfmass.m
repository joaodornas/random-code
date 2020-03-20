speed=peak_tf./peak_sf;

for i=1:6
    responsesf=eval(strcat('response_sf',int2str(i)));
    responsetf=eval(strcat('response_tf',int2str(i)));
     for k=1:105
         t=isnan(responsesf(k));
         if t==1
             responsesf(k)=0;
         end
         t=isnan(responsetf(k));
         if t==1
             responsetf(k)=0;
         end
     end
    nonzero=find(responsesf>0.00001);
    nonzero_responsesf=responsesf(nonzero);
    nonzero_peak_sf=peak_sf(nonzero);
    center_mass_sf=2*(sum(log2(nonzero_responsesf.*nonzero_peak_sf))/sum(nonzero_responsesf));
    assignin('base',strcat('center_mass_sf',int2str(i)),center_mass_sf);
    nonzero=find(responsetf>0);
    nonzero_responsetf=responsetf(nonzero);
    nonzero_peak_tf=peak_tf(nonzero);
    center_mass_tf=2*(sum(log2(nonzero_responsetf.*nonzero_peak_tf))/sum(nonzero_responsetf));
    assignin('base',strcat('center_mass_tf',int2str(i)),center_mass_tf);
    for j=1:6
       responsetf_speed=eval(strcat('response_tf',int2str(j)));
       response=responsesf.*responsetf_speed;
        for k=1:105
             t=isnan(response(k));
             if t==1
                 response(k)=0;
             end
        end
       nonzero=find(response>0);
       nonzero_response=response(nonzero);
       nonzero_speed=speed(nonzero);
       center_mass=2*(sum(log2(nonzero_response.*nonzero_speed))/sum(nonzero_response));
       assignin('base',strcat('center_mass_sf',int2str(i),'_tf',int2str(j)),center_mass);
    end
end