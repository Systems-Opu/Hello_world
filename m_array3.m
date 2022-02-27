[numerict txtt]=xlsread('geneList.xls');

%[numeric1 txt1]=xlsread('dfr1.xls');%First Set
[numeric2 txt2]=xlsread('dfr2.xls');%Set3 &set4
[numeric txt]=xlsread('dfr.xls');%knockout set


r=size(txtt);dfr1=[];dfr=[];dfr2=[];
sr=size(txt1);

for l=1:r(1)
    g=txtt(l,1);
    dfr1= [dfr1;g];
    dfr = [dfr;g];
    dfr2= [dfr2;g];
end


for l=1:r(1)
    g=txtt(l,1);
    for j=1:size(numeric1,1)%for First Set
        gf=txt1(j,1);e=[];
        if strcmp(g,gf)==1
         e=numeric1(j,2:5);
         dfr1{l,2}=e(1);
          dfr1{l,3}=e(2);
           dfr1{l,4}=e(3);
            dfr1{l,5}=e(4);
        end
    end
    
    for j=1:size(numeric,1)%for Knock Out Set
        gf=txt(j,9);e=[];
        if strcmp(g,gf)==1
         e=numeric(j,1:8);
         dfr{l,2}=e(1);
          dfr{l,3}=e(2);
           dfr{l,4}=e(3);
            dfr{l,5}=e(4);
             dfr{l,6}=e(5);
              dfr{l,7}=e(6);
               dfr{l,8}=e(7);
                dfr{l,9}=e(8);
        end
    end
    
    for j=1:size(numeric2,1)%for Set 3 &Set 4
        gf=txt2(j,1);e=[];
        if strcmp(g,gf)==1
         e=numeric2(j,2:9);
         dfr2{l,2}=e(1);
          dfr2{l,3}=e(2);
           dfr2{l,4}=e(3);
            dfr2{l,5}=e(4);
             dfr2{l,6}=e(5);
              dfr2{l,7}=e(6);
               dfr2{l,8}=e(7);
                dfr2{l,9}=e(8);
        end
    end
    
end
gd2=[];yt=[];
for t=1:size(dfr2,1)
gd2(t,1)=dfr2{t,2};
gd2(t,2)=dfr2{t,3};
gd2(t,3)=dfr2{t,4};
gd2(t,4)=dfr2{t,5};
gd2(t,5)=dfr2{t,6};
gd2(t,6)=dfr2{t,7};
gd2(t,7)=dfr2{t,8};
gd2(t,8)=dfr2{t,9};
end

gd=[];yt=[];
for t=1:size(dfr2,1)
gd(t,1)=dfr{t,2};
gd(t,2)=dfr{t,3};
gd(t,3)=dfr{t,4};
gd(t,4)=dfr{t,5};
gd(t,5)=dfr{t,6};
gd(t,6)=dfr{t,7};
gd(t,7)=dfr{t,8};
gd(t,8)=dfr{t,9};
end

gd1=[];yt=[];
for t=1:size(dfr2,1)
gd1(t,1)=dfr1{t,2};
gd1(t,2)=dfr1{t,3};
gd1(t,3)=dfr1{t,4};
gd1(t,4)=dfr1{t,5};

end

xlswrite('test data.xls',gd)
%xlswrite('test data.xls',gd1)
xlswrite('test data.xls',gd2)