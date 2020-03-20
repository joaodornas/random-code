function Und_subjects_new_parcellation
 %subjects=['CK'; 'GM'; 'GS'; 'PP'; 'TC'];
 
 folder{1} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-1-22-10-2015/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 folder{2} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-2-26-10-2015/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 folder{3} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-3-3-11-2015/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 folder{4} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-4-2-11-2015/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 folder{5} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-5-2-11-2015/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 folder{6} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-6-24-11-2015/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 folder{7} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-7-14-01-2016/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 folder{8} = '/Volumes/INDIREA/_DATA/LOW-HIGH-ATTENTION/SUBJECT-8-14-01-2016/preprocessed/B0-DTI/6.forbedpost.bedpostX-SP758/';
 
 nSubjects = 8;
 
 for iSubject=[8]
     
    disp(strcat('iSubject:',int2str(iSubject))); 
     
    subjects='xx';
    th=5;

    nNodes = 758;

    file_mat=[strcat(folder{iSubject},'TOTAL_connmatrix.txt')];
    C=load(file_mat);
    file_vec=[strcat(folder{iSubject},'ROIs_vector.txt')];
    voxvect=load(file_vec);
    voxmat=repmat(voxvect,1,nNodes);

    size(C)
    size(voxmat)

    C=C./voxmat;

    % Set diagonal to zero
    for n=1:nNodes
        C(n,n)=0;
    end

     % Create undirectional matrix
     for i=1:nNodes, 
         for j=1:nNodes, 
             new(i,j)=mean([C(i,j) C(j,i)],2); 
         end;
     end;

     C=new;

    % Threshold to th=5
    I=find(C<th);
    C(I)=0;



    file_save=[strcat(folder{iSubject},'Und_PRE_xx')];
    save(file_save,'C')

 end
 
end