

    %%% CHECK DATA ON DOCS PAPERS

    database_label = 'retina';

    conn = database(strcat('docs_papers_',database_label),'root','senha1111','Vendor','MySQL','Server','127.0.0.1','PortNumber',3306);

    %%% check PDF as Image

    iddoc = 850;
    nPages = 3;

    for iPage=1:nPages

        select = strcat('SELECT image FROM docs_papers_retina.pagesasimages WHERE iddoc = ',int2str(iddoc),' AND pageid =',int2str(iPage));

        curs = exec(conn,select);
        curs = fetch(curs);
        
        result = curs.Data';
        
        t = result{1};
        img(1,:) = t;
        
        t = result{2};
        img(2,:) = t;
        
        t = result{3};
        img(3,:) = t;
        
    end


