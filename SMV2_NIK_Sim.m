tic
%Add Model
%=========
SMV2=sbiomodel('SMV2');
%Add species 
%===========
Species={'Ikba','Ikbb','Ikbe','Ikbd','p100','p105','RelA','RelB','p50','p52','R_lB',            ...
         'RelAp50','RelAp52','RelBp50','RelBp52','RelBp100','RelAp100','R_lBp50','R_lBp100',    ...
         'IkbaRelAp50','IkbbRelAp50','IkbeRelAp50','IkbdRelAp50',                               ...
         'IkbaRelAp52','IkbbRelAp52','IkbeRelAp52','IkbdRelAp52',                               ...
         'IkbaRelBp50','IkbbRelBp50','IkbeRelBp50','IkbdRelBp50',                               ...
         'IkbaR_lBp50','IkbbR_lBp50','IkbeR_lBp50','IkbdR_lBp50',                               ...
         'tIkba','tIkbb','tIkbe','tIkbd','tp100','tp105','tRelA','tRelB','tp50','tp52','tR_lB', ...
         'NEMO_ACTIVITY','NIK_ACTIVITY'};
Speciesn=cell(1,34);
for w=1:length(Species)
if w >= 36
continue
end
S=Species{w};S(end+1)='n'; Speciesn{w}=S;
end
I='InitialAmountUnit';U='nanomole';i='InitialAmount';C='ConstantAmount';B='BoundaryCondition';t='true';f='false';v=' ';
for w=1:length(Species)
if w<=35
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
else
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),true,sprintf(i),0,sprintf(I),sprintf(U));']);
end
     if w>=36;continue;end
     eval([sprintf(Speciesn{w}) '=addspecies(SMV2,Speciesn{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
end
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'InitialAmount',1.00);
set (NEMO_ACTIVITY, 'ConstantAmount', false);
set (NEMO_ACTIVITY, 'BoundaryCondition', true);
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'ConstantAmount', false);
set (NIK_ACTIVITY, 'BoundaryCondition', true);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

%Add Reactions and Parameters
%============================
%Monomers
     %                                       1           2           3         4          5          6          7           8           9         10         11    
     Parameters=struct;
     %                                      ikba        ikbb       ikbe       ikbd       p100       p105       RelA       RelB        P50        P52      R_lB   
     Parameters.monparams.Const_trnx     = [0.43e-2     4.62e-4    8.4e-5      0          6e-4       0         1e-5       3.4e-5     1.4e-3      0         3.4e-5   ]; 
                                                                                          
     Parameters.monparams.mRNA_Deg       = [4.4e-2      6.9e-3     3.8e-3      0          1.9e-3     0         2.9e-3     2.4e-2     2.9e-3      0         2.4e-2   ];
     %                                                                                                                    
     Parameters.monparams.Trnsln_rate    = [12          12         12          0          12         0         12         12         12          0         12       ];
    
     Parameters.monparams.CnstDegFrSpCyt = [1.5e-2      2.4        1.8         1.1e-2     5.7e-3     0         5.7e-3     3.9e-3     3.7e-2      5.7e-3    3.9e-3   ];
     %                                                                                                                    
     Parameters.monparams.CnstDegFrSpNuc = [1.5e-2      2.4        1.8         1.1e-2     5.7e-3     0         5.7e-3     3.9e-3     3.7e-2      5.7e-3    3.9e-3   ];
     %                                                                                                                    
     Parameters.monparams.NEMODegFrSp    = [1.4e-3      4.5e-4     9e-4        0          0          0         0          0          0           0         0        ];
      
     Parameters.monparams.NIKDegFrSp     = [0           0          0           1.2        0          0         0          0          0           0         0        ];
     
     Parameters.monparams.Nuc_Imp        = [6e-2        9e-3       4.5e-2      4.5e-2     0          0         0          0          0           0         0        ];
	
     Parameters.monparams.Nuc_Exp        = [1.2e-2      1.2e-2     1.2e-2      1.2e-2     0          0         0          0          0           0         0        ];    
  %========================================================================
     
     %Dimers                                1       2        3          4         5         6         7           8
     %                                    A50      A52      B50        B52       BP100     AP100     R_lBp50    R_lBp100
     Parameters.dimparams.MonAssCyt    = [9.9e-4   9.6e-4   9.9e-4     9.6e-4    1.7e-2    5.2e-4     9.9e-4      1.7e-2  ]; 
     %                                  
     Parameters.dimparams.MonAssNuc    = [9.9e-4   9.6e-4   9.9e-4     9.6e-4    0         0          9.9e-4      0       ];
     %                                  
     Parameters.dimparams.DimDissCyt   = [1.4e-3   1.4e-3   1.8e-3     1.4e-3    1.1e-3    1.1e-3     1.8e-3      1.1e-3  ];
     %                                                                 4.4
     Parameters.dimparams.DimDissNuc   = [1.4e-3   1.4e-3   1.8e-3     1.4e-3    0         0          1.8e-3      0       ];
     %                                                                 4.4
     Parameters.dimparams.NucImp       = [5.4      5.4      5.4        5.4       4.8e-3    4.8e-3     5.4         4.8e-3  ];
     
     Parameters.dimparams.NucExp       = [4.8e-3   4.8e-3   4.8e-3     4.8e-3    4.8e-3    4.8e-3     4.8e-3      4.8e-3  ];
     
     Parameters.dimparams.DimDegCyt    = [2.4e-4   2.4e-4   2.3e-4     4.6e-4    5e-3      2.4e-4     2.3e-4      5e-3    ];
     %                                                                 9.6
     Parameters.dimparams.DimDegNuc    = [2.4e-4   2.4e-4   2.3e-4     4.6e-4    5e-3      2.4e-4     2.3e-4      5e-3    ];
     %                                                                 9.6
%======================================   
%===================================================  
%Trimer                                                 1         2         3          4         5         6          7          8           9           10         11          12        13         14         15         16
%                                                     RelA50:a  RelA50:b  RelA50:e   RelA50:d  RelA52:a  RelA52:b  RelA52:e   RelA52:d    RelB50:a    RelB50:b    RelB50:e    RelB50:d   RelB_50:a  RelB_50:b  RelB_50:e  RelB_50:d
Parameters.trimparams.Trim_Ass_Cyt                   = [3e-1      2e-1     2e-1       2e-2      2e-1      2e-1      2e-1       2e-2        9e-2          0         4e-2        1.2e-1     0.001        0      0.001       1.2e-1  ];   % Association for composite species in Cytoplasm, NFkB + IkB -> NFkB:IkB
%                                                                                                                                                              
Parameters.trimparams.Trim_Ass_Nuc                   = [3e-1      2e-1     2e-1       2e-2      2e-1      2e-1      2e-1       2e-2        9e-2          0         4e-2        1.2e-1     0.001        0      0.001       1.2e-1  ];   % Association for composite species in Nucleus, NFkBn + IkBn -> NFkB:IkBn
%                                                                                                                          
Parameters.trimparams.Trim_Dis_Cyt                   = [8.4e-3    3.4e-2   8.4e-3     8.4e-2    8.4e-3    3.4e-2    8.4e-3     8.4e-2      8.4e-3        0         8.4e-3      1.4e-3     8.4e-2       0      8.4e-2      1.4e-3  ];   % Dissociation for composite species, NFkB:IkB -> NFkB + IkB
%                                                                                                                                       
Parameters.trimparams.Trim_Dis_Nuc                   = [8.4e-3    3.4e-2   8.4e-3     8.4e-2    8.4e-3    3.4e-2    8.4e-3     8.4e-2      8.4e-3        0         8.4e-3      1.4e-3     8.4e-2       0      8.4e-2      1.4e-3  ];   % Dissociation for composite species, NFkB:IkB -> NFkB + IkB  
%                                                                                                                                      
Parameters.trimparams.Const_Trim_Deg_IkB_Cyt         = [2e-4      2e-4     2e-4       2e-4      2e-4      2e-4      2e-4       2e-4        2e-4          0         2e-4        2e-4       2e-4         0      2e-4        2e-4    ];   % Const. deg IkBs within composite species NFkB:IkB -> NFkB

Parameters.trimparams.Const_Trim_Deg_IkB_Nuc         = [2e-4      2e-4     2e-4       2e-4      2e-4      2e-4      2e-4       2e-4        2e-4          0         2e-4        2e-4       2e-4         0      2e-4        2e-4    ];   % Const. deg IkBs within composite species NFkB:IkB -> NFkB
 
Parameters.trimparams.Const_Trim_Deg_NFkb_Cyt        = [2e-4      2e-4     2e-4       2e-4      2e-4      2e-4      2e-4       2e-4        3.6e-4        0         2.4e-4      2.4e-4     3.6e-4       0      2.4e-4      2.4e-4  ];   % Const. deg NFkBs within composite species NFkB:IkB -> IkB

Parameters.trimparams.Const_Trim_Deg_NFkb_Nuc        = [2e-4      2e-4     2e-4       2e-4      2e-4      2e-4      2e-4       2e-4        3.6e-4        0         2.4e-4      2.4e-4     3.6e-4       0      2.4e-4      2.4e-4  ];   % Const. deg NFkBs within composite species NFkB:IkB -> IkB

Parameters.trimparams.NEMO_Deg_Ikb_inTrim            = [1.4e-3    4.5e-4   9e-4       0         1.4e-3    4.5e-4    9e-4       0           1.4e-3        0         9e-4        0          1.4e-3       0      9e-4        0       ];   % NEMO mediated deg of IkBs within composite species
%                                                                                                                          
Parameters.trimparams.NIK_Deg_Ikb_inTrim             = [0         0        0          1e-3      0         0         0          1e-3        0             0         0           1e-3       0            0      0           1e-3    ];   % NIK mediated deg of IkBs within composite species

Parameters.trimparams.Nuc_Imp_Trim                   = [2.8e-1    2.8e-2   1.4e-1     2.8e-2    2.8e-1    2.8e-2    1.4e-1     2.8e-2      2.8e-1        0         1.4e-1      2.8e-2     2.8e-1       0      1.4e-1      2.8e-2  ];   % Import Reactions of IkB's:ReLAP50/ReLAP52

Parameters.trimparams.Nuc_Exp_Trim                   = [8.4e-1    4.2e-1   4.2e-1     2.8e-1    8.4e-1    4.2e-1    4.2e-1     4.2e-1      8.4e-1        0         4.2e-1      4.2e-1     8.4e-1       0      4.2e-1      4.2e-1  ];   % Export Reactions of IkB's:ReLAP50/ReLAP52

%====================================================
                    %RelAp50 RelAp52 RelBp50  R_lBp50              ||   Mirrored Boxes:
        Parameters.W = [20     20      0        0     ;... %Ikba   ||   W11      W12      W13      W14     25     25      0    0
%                                                                  ||
                        15     15      0        0     ;... %Ikbb   ||   W21      W22      W23      W24     10     10      0    0
%                                                                  ||
                        25     25      0        0     ;... %Ikbe   ||   W31      W32      W33      W34     125    125     0    0
%                                                                  ||                                          
                        30     30      0        0     ;... %P100   ||   W41      W42      W43      W44     50     50      0    0
%                                                                  ||                        
                        9      9       9        9     ;... %RelB   ||   W51      W52      W53      W54     30     30      50   50
%                                                                  ||
                        9      9       9        9  ]  ;    %R_lB   ||   W61      W62      W63      W64     30     30      50   50
%                                                                  ||
%                     =============================================||==================================    ==================
       Parameters.Kd = [ 85     85     0       0      ;... %Ikba   ||   Kd11     Kd12     Kd13     Kd14     125    125     0    0     
%                                                                  ||
                         85     85     0       0      ;... %Ikbb   ||   Kd21     Kd22     Kd23     Kd24     150    150     0    0    
%                                                                  ||
                         85     85     0       0      ;... %Ikbe   ||   Kd31     Kd32     Kd33     Kd34     150    150     0    0       
%                                                                  ||                                    
                         85     85     0       0      ;... %P100   ||   Kd41     Kd42     Kd43     Kd44     100    100     0    0    
%                                                                  ||                       
                         10     10     3.33    3.33     ;... %RelB ||   Kd51     Kd52     Kd53     Kd54     170    170     45   45
%                                                                  ||
                         10     10     3.33    3.33  ]  ;    %R_lB ||   Kd61     Kd62     Kd63     Kd64     170    170     45   45
%=========================================================
%Special Reactions
   Parameters.p100_Ass            = 9e-2;     % p100 association rate to form IkBd in cytoplasm %%% 1.2E-3
   %                           K253 
   Parameters.Ikbd_Diss           = 1.2e-5;   % IkBd dissocation rate to form p100 in cytoplasm %%% 1.2E-2
   %                           K254
   Parameters.NIK_mediated_p52    = 3.5e-3;   % NIK mediated generation of P52 through processing
   %                           K255
%                       
M='MassAction';n=1;V='ValueUnits'; U0='nanomole/minute';U1='1/minute';U2='(1/nanomole)*(1/minute)';U3='1/nanomole';
for q1=1:11
  if (q1~=8) && (q1~=11)
      %Constitutive transcription
      Rkonst=horzcat('Kct',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst),Parameters.monparams.Const_trnx(q1),V,U0);']);n=n+1; 
  
      %Contitutive mRNA degradation
      Rconst=horzcat('Kmdeg',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.monparams.mRNA_Deg(q1),V,U1);']);n=n+1;
  end
  if (q1~=4) && (q1~=10)
      %Translation Rate
      R3=horzcat(Species{q1+35},' ','->',' ',Species{q1});RO3=horzcat('Trnsln_',Species{q1},num2str(n));KLO3=horzcat('kine',num2str(n));Rkonst3=horzcat('K',num2str(n));
      eval([sprintf(RO3) '=addreaction(SMV2,sprintf(R3));']); 
      eval([sprintf(KLO3) '= addkineticlaw(eval(sprintf(RO3)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst3), Parameters.monparams.Trnsln_rate (q1),V,U1);']);
      set(eval(sprintf(KLO3)),'ParameterVariableName', {sprintf(Rkonst3)});n=n+1;
  end
      %Constitutive Degradation of Free primary Species (Cytoplasm)
      R4=horzcat(Species{q1},' ','->',' ','null');RO4=horzcat('Const_Prot_Deg_',Species{q1},num2str(n));KLO4=horzcat('kine',num2str(n));Rkonst4=horzcat('K',num2str(n));
      eval([sprintf(RO4) '=addreaction(SMV2,sprintf(R4));']); 
      eval([sprintf(KLO4) '= addkineticlaw(eval(sprintf(RO4)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst4),Parameters.monparams.CnstDegFrSpCyt(q1),V,U1);']); 
      set(eval(sprintf(KLO4)),'ParameterVariableName', {sprintf(Rkonst4)});n=n+1;
      %Constitutive Degradation of Free primary Species (Nucleus)
      R5=horzcat(Speciesn{q1},' ','->',' ','null');RO5=horzcat('Const_Prot_Deg_',Speciesn{q1},num2str(n));KLO5=horzcat('kine',num2str(n));Rkonst5=horzcat('K',num2str(n));
      eval([sprintf(RO5) '=addreaction(SMV2,sprintf(R5));']); 
      eval([sprintf(KLO5) '= addkineticlaw(eval(sprintf(RO5)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst5),Parameters.monparams.CnstDegFrSpNuc(q1),V,U1);']); 
      set(eval(sprintf(KLO5)),'ParameterVariableName', {sprintf(Rkonst5)});n=n+1;
      %NEMO Mediated Degradation of Canonical Ikbs
      if q1<=3
      R6=horzcat(Species{q1},' ','+',' ','NEMO_ACTIVITY',' ','->',' ','null');RO6=horzcat('NEMO_Mon_Deg_',Species{q1},num2str(n));KLO6=horzcat('kine',num2str(n));Rkonst6=horzcat('K',num2str(n));
      eval([sprintf(RO6) '=addreaction(SMV2,sprintf(R6));']); 
      eval([sprintf(KLO6) '= addkineticlaw(eval(sprintf(RO6)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst6),Parameters.monparams.NEMODegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO6)),'ParameterVariableName', {sprintf(Rkonst6)});n=n+1;
      end
      %NIK Mediated Degradation of Non-Canonical Ikb
      if q1==4
      R7=horzcat(Species{q1},' ','+',' ','NIK_ACTIVITY',' ','->',' ','null');RO7=horzcat('NIK_Mon_Deg_',Species{q1},num2str(n));KLO7=horzcat('kine',num2str(n));Rkonst7=horzcat('K',num2str(n));
      eval([sprintf(RO7) '=addreaction(SMV2,sprintf(R7));']); 
      eval([sprintf(KLO7) '= addkineticlaw(eval(sprintf(RO7)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst7), Parameters.monparams.NIKDegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO7)),'ParameterVariableName', {sprintf(Rkonst7)});n=n+1;
      end
      %Nuclear Import 
      if q1<=4
      R8=horzcat(Species{q1},' ','->',' ',Speciesn{q1});RO8=horzcat('Nuc_Im_Mon',Species{q1},num2str(n));KLO8=horzcat('kine',num2str(n));Rkonst8=horzcat('K',num2str(n));
      eval([sprintf(RO8) '=addreaction(SMV2,sprintf(R8));']); 
      eval([sprintf(KLO8) '= addkineticlaw(eval(sprintf(RO8)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst8),Parameters.monparams.Nuc_Imp(q1),V,U1);']); 
      set(eval(sprintf(KLO8)),'ParameterVariableName', {sprintf(Rkonst8)});n=n+1;
      end
      %Nuclear Export
      if q1<=4
      R9=horzcat(Speciesn{q1},' ','->',' ',Species{q1});RO9=horzcat('Nuc_Ex_Mon',Species{q1},num2str(n));KLO9=horzcat('kine',num2str(n));Rkonst9=horzcat('K',num2str(n));
      eval([sprintf(RO9) '=addreaction(SMV2,sprintf(R9));']); 
      eval([sprintf(KLO9) '= addkineticlaw(eval(sprintf(RO9)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst9),Parameters.monparams.Nuc_Exp(q1),V,U1);']); 
      set(eval(sprintf(KLO9)),'ParameterVariableName', {sprintf(Rkonst9)});n=n+1;
      end
end
%==================================================== 


for q2=1:8
       %Association  of Dimers in Cytoplasm
       R10=horzcat(Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end),' ','->',' ',Species{11+q2});RO10=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO10=horzcat('kine',num2str(n));Rkonst10=horzcat('K',num2str(n));
       eval([sprintf(RO10) '=addreaction(SMV2,sprintf(R10));']); 
       eval([sprintf(KLO10) '= addkineticlaw(eval(sprintf(RO10)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst10),Parameters.dimparams.MonAssCyt(q2),V,U2);']); 
       set(eval(sprintf(KLO10)),'ParameterVariableName', {sprintf(Rkonst10)});n=n+1;
       %Association  of Dimers in Nucleus
       R11=horzcat(Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end),' ','->',' ',Speciesn{11+q2});RO11=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO11=horzcat('kine',num2str(n));Rkonst11=horzcat('K',num2str(n));
       eval([sprintf(RO11) '=addreaction(SMV2,sprintf(R11));']); 
       eval([sprintf(KLO11) '= addkineticlaw(eval(sprintf(RO11)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst11),Parameters.dimparams.MonAssNuc(q2),V,U2);']); 
       set(eval(sprintf(KLO11)),'ParameterVariableName', {sprintf(Rkonst11)});n=n+1;
       %Dissociation  of Dimers in Cytoplasm
       R12=horzcat(Species{11+q2},' ','->',' ',Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end));RO12=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO12=horzcat('kine',num2str(n));Rkonst12=horzcat('K',num2str(n));
       eval([sprintf(RO12) '=addreaction(SMV2,sprintf(R12));']); 
       eval([sprintf(KLO12) '= addkineticlaw(eval(sprintf(RO12)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst12),Parameters.dimparams.DimDissCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO12)),'ParameterVariableName', {sprintf(Rkonst12)});n=n+1;
       %Dissociation  of Dimers in Nucleus
       R13=horzcat(Speciesn{11+q2},' ','->',' ',Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end));RO13=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO13=horzcat('kine',num2str(n));Rkonst13=horzcat('K',num2str(n));
       eval([sprintf(RO13) '=addreaction(SMV2,sprintf(R13));']); 
       eval([sprintf(KLO13) '= addkineticlaw(eval(sprintf(RO13)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst13),Parameters.dimparams.DimDissNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO13)),'ParameterVariableName', {sprintf(Rkonst13)});n=n+1;
       %Nuclear Import
       R14=horzcat(Species{11+q2},' ','->',' ',Speciesn{11+q2});RO14=horzcat('Nuc_Im_Dim',Species{q2+11},num2str(n));KLO14=horzcat('kine',num2str(n));Rkonst14=horzcat('K',num2str(n));
       eval([sprintf(RO14) '=addreaction(SMV2,sprintf(R14));']); 
       eval([sprintf(KLO14) '= addkineticlaw(eval(sprintf(RO14)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst14),Parameters.dimparams.NucImp(q2),V,U1);']); 
       set(eval(sprintf(KLO14)),'ParameterVariableName', {sprintf(Rkonst14)});n=n+1;
       %Nuclear Export
       R15=horzcat(Speciesn{11+q2},' ','->',' ',Species{11+q2});RO15=horzcat('Nuc_Ex_Dim',Species{q2+11},num2str(n));KLO15=horzcat('kine',num2str(n));Rkonst15=horzcat('K',num2str(n));
       eval([sprintf(RO15) '=addreaction(SMV2,sprintf(R15));']); 
       eval([sprintf(KLO15) '= addkineticlaw(eval(sprintf(RO15)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst15),Parameters.dimparams.NucExp(q2),V,U1);']); 
       set(eval(sprintf(KLO15)),'ParameterVariableName', {sprintf(Rkonst15)});n=n+1;
       %Free Dimer Degradation in Cytoplasm
       R16=horzcat(Species{11+q2},' ','->',' ','null');RO16=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO16=horzcat('kine',num2str(n));Rkonst16=horzcat('K',num2str(n));
       eval([sprintf(RO16) '=addreaction(SMV2,sprintf(R16));']); 
       eval([sprintf(KLO16) '= addkineticlaw(eval(sprintf(RO16)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst16),Parameters.dimparams.DimDegCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO16)),'ParameterVariableName', {sprintf(Rkonst16)});n=n+1;
       %Free Dimer Degradation in Nucleus
       R17=horzcat(Speciesn{11+q2},' ','->',' ','null');RO17=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO17=horzcat('kine',num2str(n));Rkonst17=horzcat('K',num2str(n));
       eval([sprintf(RO17) '=addreaction(SMV2,sprintf(R17));']); 
       eval([sprintf(KLO17) '= addkineticlaw(eval(sprintf(RO17)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst17),Parameters.dimparams.DimDegNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO17)),'ParameterVariableName', {sprintf(Rkonst17)});n=n+1;
end
%========================================================
                                                     
for q3=1:16
        %Association of trimers in Cytoplasm
        R18=horzcat(Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end),' ','->',' ',Species{q3+19});RO18=horzcat('Trim_Ass_Cyto',Species{q3+19},num2str(n));KLO18=horzcat('kine',num2str(n));Rkonst18=horzcat('K',num2str(n));
        eval([sprintf(RO18) '=addreaction(SMV2,sprintf(R18));']); 
        eval([sprintf(KLO18) '= addkineticlaw(eval(sprintf(RO18)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst18),Parameters.trimparams.Trim_Ass_Cyt(q3),V,U2);']); 
        set(eval(sprintf(KLO18)),'ParameterVariableName', {sprintf(Rkonst18)});n=n+1;
        %Association of trimers in Nucleus
        R19=horzcat(Speciesn{q3+19}(1:4),'n',' ','+',' ',Speciesn{q3+19}(5:end),' ','->',' ',Speciesn{q3+19});RO19=horzcat('Trim_Ass_Nuc',Species{q3+19},num2str(n));KLO19=horzcat('kine',num2str(n));Rkonst19=horzcat('K',num2str(n));
        eval([sprintf(RO19) '=addreaction(SMV2,sprintf(R19));']); 
        eval([sprintf(KLO19) '= addkineticlaw(eval(sprintf(RO19)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst19),Parameters.trimparams.Trim_Ass_Nuc(q3),V,U2);']); 
        set(eval(sprintf(KLO19)),'ParameterVariableName', {sprintf(Rkonst19)});n=n+1;
        %Dissociation of trimers in Cytoplasm
        R20=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end));RO20=horzcat('Trim_Dis_Cyt',Species{q3+19},num2str(n));KLO20=horzcat('kine',num2str(n));Rkonst20=horzcat('K',num2str(n));
        eval([sprintf(RO20) '=addreaction(SMV2,sprintf(R20));']); 
        eval([sprintf(KLO20) '= addkineticlaw(eval(sprintf(RO20)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst20),Parameters.trimparams.Trim_Dis_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO20)),'ParameterVariableName', {sprintf(Rkonst20)});n=n+1;
        %Dissociation of trimers in Nucleus
        R21=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n',' ','+',' ',Species{q3+19}(5:end),'n');RO21=horzcat('Trim_Dis_Nuc',Species{q3+19},num2str(n));KLO21=horzcat('kine',num2str(n));Rkonst21=horzcat('K',num2str(n));
        eval([sprintf(RO21) '=addreaction(SMV2,sprintf(R21));']); 
        eval([sprintf(KLO21) '= addkineticlaw(eval(sprintf(RO21)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst21),Parameters.trimparams.Trim_Dis_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO21)),'ParameterVariableName', {sprintf(Rkonst21)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Cytoplasm
        R22=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(5:end));RO22=horzcat('Const_Deg_Trim_Ikb_Cyt',Species{q3+19},num2str(n));KLO22=horzcat('kine',num2str(n));Rkonst22=horzcat('K',num2str(n));
        eval([sprintf(RO22) '=addreaction(SMV2,sprintf(R22));']); 
        eval([sprintf(KLO22) '= addkineticlaw(eval(sprintf(RO22)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst22),Parameters.trimparams.Const_Trim_Deg_IkB_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO22)),'ParameterVariableName', {sprintf(Rkonst22)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Nucleus
        R23=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(5:end),'n');RO23=horzcat('Const_Deg_Trim_Ikb_Nuc',Species{q3+19},num2str(n));KLO23=horzcat('kine',num2str(n));Rkonst23=horzcat('K',num2str(n));
        eval([sprintf(RO23) '=addreaction(SMV2,sprintf(R23));']); 
        eval([sprintf(KLO23) '= addkineticlaw(eval(sprintf(RO23)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst23),Parameters.trimparams.Const_Trim_Deg_IkB_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO23)),'ParameterVariableName', {sprintf(Rkonst23)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Cytoplasm
        R24=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4));RO24=horzcat('Const_Deg_Trim_NFkb_Cyt',Species{q3+19},num2str(n));KLO24=horzcat('kine',num2str(n));Rkonst24=horzcat('K',num2str(n));
        eval([sprintf(RO24) '=addreaction(SMV2,sprintf(R24));']); 
        eval([sprintf(KLO24) '= addkineticlaw(eval(sprintf(RO24)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst24),Parameters.trimparams.Const_Trim_Deg_NFkb_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO24)),'ParameterVariableName', {sprintf(Rkonst24)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Nucleus
        R25=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n');RO25=horzcat('Const_Deg_Trim_NFkb_Nuc',Species{q3+19},num2str(n));KLO25=horzcat('kine',num2str(n));Rkonst25=horzcat('K',num2str(n));
        eval([sprintf(RO25) '=addreaction(SMV2,sprintf(R25));']); 
        eval([sprintf(KLO25) '= addkineticlaw(eval(sprintf(RO25)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst25),Parameters.trimparams.Const_Trim_Deg_NFkb_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO25)),'ParameterVariableName', {sprintf(Rkonst25)});n=n+1;
        %NEMO Mediated Degradation of Ikb in Trimer in Cytoplasm
        R26=horzcat(Species{q3+19},' ','+',' ','NEMO_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO26=horzcat('NEMO_DEg_inTrim',Species{q3+19},num2str(n));KLO26=horzcat('kine',num2str(n));Rkonst26=horzcat('K',num2str(n));
        eval([sprintf(RO26) '=addreaction(SMV2,sprintf(R26));']); 
        eval([sprintf(KLO26) '= addkineticlaw(eval(sprintf(RO26)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst26),Parameters.trimparams.NEMO_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO26)),'ParameterVariableName', {sprintf(Rkonst26)});n=n+1;
        %NIK Mediated Degradation of Ikb in Trimer in Cytoplasm
        R27=horzcat(Species{q3+19},' ','+',' ','NIK_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO27=horzcat('NIK_DEg_inTrim',Species{q3+19},num2str(n));KLO27=horzcat('kine',num2str(n));Rkonst27=horzcat('K',num2str(n));
        eval([sprintf(RO27) '=addreaction(SMV2,sprintf(R27));']); 
        eval([sprintf(KLO27) '= addkineticlaw(eval(sprintf(RO27)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst27),Parameters.trimparams.NIK_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO27)),'ParameterVariableName', {sprintf(Rkonst27)});n=n+1;
        %Nuclear Import of Trimers
        R28=horzcat(Species{q3+18},' ','->',Speciesn{q3+18});RO28=horzcat('Nuc_Imp_Trim',Species{q3+18},num2str(n));KLO28=horzcat('kine',num2str(n));Rkonst28=horzcat('K',num2str(n));
        eval([sprintf(RO28) '=addreaction(SMV2,sprintf(R28));']); 
        eval([sprintf(KLO28) '= addkineticlaw(eval(sprintf(RO28)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst28),Parameters.trimparams.Nuc_Imp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO28)),'ParameterVariableName', {sprintf(Rkonst28)});n=n+1;
        %Nuclear Export of Trimers
        R29=horzcat(Speciesn{q3+19},' ','->',Species{q3+19});RO29=horzcat('Nuc_Exp_Trip',Species{q3+19},num2str(n));KLO29=horzcat('kine',num2str(n));Rkonst29=horzcat('K',num2str(n));
        eval([sprintf(RO29) '=addreaction(SMV2,sprintf(R29));']); 
        eval([sprintf(KLO29) '= addkineticlaw(eval(sprintf(RO29)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst29),Parameters.trimparams.Nuc_Exp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO29)),'ParameterVariableName', {sprintf(Rkonst29)});n=n+1;
end

   
        R30=horzcat('2',' ',Species{5},' ','->',Species{4});RO30=horzcat('p100_Ass',Species{5},num2str(n));KLO30=horzcat('kine',num2str(n));Rkonst30=horzcat('K',num2str(n));
        eval([sprintf(RO30) '=addreaction(SMV2,sprintf(R30));']); 
        eval([sprintf(KLO30) '= addkineticlaw(eval(sprintf(RO30)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst30),Parameters.p100_Ass,V,U2);']); 
        set(eval(sprintf(KLO30)),'ParameterVariableName', {sprintf(Rkonst30)});n=n+1;

        R31=horzcat(Species{4},' ','->',' ','2',' ',Species{5});RO31=horzcat('p100_Ass',Species{4},num2str(n));KLO31=horzcat('kine',num2str(n));Rkonst31=horzcat('K',num2str(n));
        eval([sprintf(RO31) '=addreaction(SMV2,sprintf(R31));']); 
        eval([sprintf(KLO31) '= addkineticlaw(eval(sprintf(RO31)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst31),Parameters.Ikbd_Diss,V,U1);']); 
        set(eval(sprintf(KLO31)),'ParameterVariableName', {sprintf(Rkonst31)});n=n+1;

        R32=horzcat(Species{5},' ','+',' ',Species{end},' ','->',' ',Species{10});RO32=horzcat('p52_processing',Species{5},num2str(n));KLO32=horzcat('kine',num2str(n));Rkonst32=horzcat('K',num2str(n));
        eval([sprintf(RO32) '=addreaction(SMV2,sprintf(R32));']); 
        eval([sprintf(KLO32) '= addkineticlaw(eval(sprintf(RO32)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst32),Parameters.NIK_mediated_p52,V,U2);']); 
        set(eval(sprintf(KLO32)),'ParameterVariableName', {sprintf(Rkonst32)});n=n+1;
%===========================================================
%Chequer Boxes of Inducibility 
%=============================
%                   
%Add Rules
%==========

RfntIkba = 'tIkba = ( Kct1 * ( ( 1 + ( W11 * ( RelAp50n / Kd11 ) ^ 3 ) + ( W12 * ( RelAp52n / Kd12 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd11 ) ^ 3 + ( RelAp52n / Kd12 ) ^ 3 ) ) ) - Kmdeg1 * tIkba';
              
RfntIkbb = 'tIkbb = ( Kct2 * ( ( 1 + ( W21 * ( RelAp50n / Kd21 ) ^ 3 ) + ( W22 * ( RelAp52n / Kd22 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd21 ) ^ 3 + ( RelAp52n / Kd22 ) ^ 3 ) ) ^ Kdel1 ) - Kmdeg3 * tIkbb';
                 
RfntIkbe = 'tIkbe = ( Kct3 * ( ( 1 + ( W31 * ( RelAp50n / Kd31 ) ^ 3 ) + ( W32 * ( RelAp52n / Kd32 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd31 ) ^ 3 + ( RelAp52n / Kd32 ) ^ 3 ) ) ^ Kdel2 ) - Kmdeg3 * tIkbe';
                
Rfntp100 = 'tp100 = ( Kct5 * ( ( 1 + ( W41 * ( RelAp50n / Kd41 ) ^ 3 ) + ( W42 * ( RelAp52n / Kd42 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd41 ) ^ 3 + ( RelAp52n / Kd42 ) ^ 3 ) ) ^ Kdel3 ) - Kmdeg5 * tp100';
                 
RfntRelA = 'tRelA = Kct7  - Kmdeg7 * tRelA';

RfntRelB = 'tRelB = ( Kct8 * ( ( 1 + ( W51 * ( RelAp50n * Kdel5 / Kd51 ) ^ 3 ) + ( W52 * ( RelAp52n * Kdel5 / Kd52 ) ^ 3 ) + ( W53 * ( RelBp50n / Kd53 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd51 ) ^ 3 + ( RelAp52n * Kdel5 / Kd52 ) ^ 3  + ( RelBp50n / Kd53 ) ^ 3 ) ) ) - Kmdeg8 * tRelB';

%RfntRelB = 'tRelB = ( Kct8 * ( ( 1 + ( W51 * ( RelAp50n / Kd51 ) ^ 3 ) + ( W52 * ( RelAp52n / Kd52 ) ^ 3 ) + ( W53 * ( RelBp50n / Kd53 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd51 ) ^ 3 + ( RelAp52n / Kd52 ) ^ 3  + ( RelBp50n / Kd53 ) ^ 3 ) ) ^ Kdel4 ) - Kmdeg8 * tRelB';

RfntR_lB = 'tR_lB = ( Kct11 * ( ( 1 + ( W61 * ( RelAp50n * Kdel5 / Kd61 ) ^ 3 ) + ( W62 * ( RelAp52n * Kdel5 / Kd62 ) ^ 3 ) + ( W63 * ( RelBp50n / Kd63 ) ^ 3 ) + ( W64 * ( R_lBp50n / Kd64 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd61 ) ^ 3 + ( RelAp52n * Kdel5 / Kd62 ) ^ 3  + ( RelBp50n / Kd63 ) ^ 3 + ( R_lBp50n / Kd64 ) ^ 3 ) ) ^ Kdel4 ) - Kmdeg11 * tR_lB';
                   
Rfntp50  = 'tp50  = Kct9 - Kmdeg9 * tp50';
for krow= 1:size(Parameters.Kd,1)
  for kcolumn= 1:size(Parameters.Kd,2)
      Rconst=horzcat('Kd',num2str(krow),num2str(kcolumn)); Rconst1=horzcat('W',num2str(krow),num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.Kd(krow,kcolumn),V,U);']);n=n+1;
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst1),Parameters.W(krow,kcolumn));']);n=n+1; 
  end
end
%============================================
  RlO1 = addrule(SMV2,sprintf(RfntIkba),'rate');
  RlO2 = addrule(SMV2,sprintf(RfntIkbb),'rate');
  RlO3 = addrule(SMV2,sprintf(RfntIkbe),'rate');
  RlO4 = addrule(SMV2,sprintf(Rfntp100),'rate');
  RlO5 = addrule(SMV2,sprintf(RfntRelA),'rate');
  RlO6 = addrule(SMV2,sprintf(RfntRelB),'rate');
  RlO7 = addrule(SMV2,sprintf(RfntR_lB),'rate');
  RlO8 = addrule(SMV2,sprintf(Rfntp50),'rate');
%============================================
Ca='ConstantValue';%VA='Value'

%DELAY
%=====
%kdels are 1 during euilibration
for kcolumn= 1:5
      Rconst=horzcat('Kdel',num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),1,Ca,false);']);
      n=n+1; 
end
%kdels are 0 at stimulation 0 min.
for P=1:4
 E=horzcat('Eee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(0));
                                         T=horzcat('time >= ',num2str(10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%kdels are 1 at stimulation phase once delay is over.
DELAY=[37 37 45 0]; %beta epsilon p100 RelB 
for P=1:4
 E=horzcat('Eeee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(1));
                                         T=horzcat('time >= ',num2str(DELAY(P)+10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%RelB star issue
%===============

      RconsT=horzcat('Kct',num2str(8));RconsTd=horzcat('Kmdeg',num2str(8));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsT),Parameters.monparams.Const_trnx(8),Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsTd),Parameters.monparams.mRNA_Deg(8),Ca,false);']); n=n+1; 
      EEE1=addevent(SMV2,'time >= 10000','Kct8 = 0');%EEE2=addevent(SMV2,'time >= 10000','Kmdeg8 = 0');

      RconsTT=horzcat('Kct',num2str(11));RconsTTd=horzcat('Kmdeg',num2str(11));
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTT),0,Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTTd),0,Ca,false);']); n=n+1; 
      %EEE3=addevent(SMV2,'time >= 150000','Kct11 =  3.6e-5');EEE4=addevent(SMV2,'time >= 150000','Kmdeg11 = 4.6e-3');
      E=horzcat('EEE',num2str(3));N=horzcat('Kct',num2str(11),' ','=',' ',num2str(Parameters.monparams.Const_trnx(11))); T=horzcat('time >= ',num2str(10000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);       
      E=horzcat('EEE',num2str(4));N=horzcat('Kmdeg',num2str(11),' ','=',' ',num2str(Parameters.monparams.mRNA_Deg(11))); %T=horzcat('time >= ',num2str(150000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);                                                                                                

%Effective RelA heterodimer concentration in nucleus involved with RelB synthesis
evalues =  [1 .01];etime = [0 60];
 for Q=0:1:60
 Kdel5 = interp1(etime, evalues, Q,'pchip');E=horzcat('Eeff',num2str(Q));N=horzcat('Kdel5',' ','=',' ',num2str(Kdel5));
                                         T=horzcat('time >= ',num2str(10030+Q));%RelA effective concentration start going down after 30 mins
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end     
      
%SIGNAL
%======
%LTb
values = [1.0000 1.4362 2.4395 2.6754 11.5747 11.4070 9.4359 9.5987];
time   = [0      30     60     180    300     480     960    1440];
%======================================
%Add Parameters of Inducibility
%==============================
for Q=0:1:1440
 NIK = interp1(time, values, Q,'pchip');E=horzcat('E',num2str(Q));N=horzcat('NIK_ACTIVITY',' ','=',' ',num2str(NIK));
                                         T=horzcat('time >= ',num2str(10000+Q));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end

%SIMULATE MODEL
%==============
CnfgstO = addconfigset(SMV2, 'De3');
set(CnfgstO.CompileOptions, 'UnitConversion', false);
set(CnfgstO.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO,'TimeUnits',' ');
set(CnfgstO,'SolverType','sundials');
set(CnfgstO,'StopTime', 11440);
sbioaccelerate(SMV2,CnfgstO)
[T1,X1,names1]= sbiosimulate(SMV2,CnfgstO);
%=========================================================
Time1 = 10000;
Time2 = 11440;

%Free Absolute Nuclear NFkB
subplot(2,1,1)
plot(T1,X1(:,24)+X1(:,26),'r','linewidth',2);hold on   %RelAp50n
plot(T1,X1(:,28)+X1(:,36)+X1(:,30),'g','linewidth',3);hold on  %RelBp50n+R_lBp50n+RelB52n   
plot(T1,X1(:,28),'c','linewidth',1);hold on %RelBp50n 
plot(T1,X1(:,36),'b','linewidth',1);hold on %R_lBp50n
plot(T1,X1(:,30),'k','linewidth',1);hold off %R_lBp50n

xlim([Time1 Time2])
ylim([0 150])
set(gca,'xtick',[0:120:1440]+10000)
%set(gca,'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' '])
set(gca,'XGrid','on')

%Transcript fold induction
for D=1:size(T)
    if T(D)>= 10000
       break
    end
end

subplot(2,1,2)

plot(T1,X1(:,71)/X1(D,71),'b');hold on   %tIkba
plot(T1,X1(:,72)/X1(D,72),'c');hold on   %tIkbb
plot(T1,X1(:,73)/X1(D,73),'m');hold on   %tIkbe
plot(T1,X1(:,75)/X1(D,75),'k');hold on   %tp100
plot(T1,X1(:,77)/X1(D,77),'r');hold on   %tRelA
plot(T1,((X1(:,78)+X1(:,81))/(X1(D,78)+X1(D,81))),'g');hold on   %tRelB
plot(T1,X1(:,79)/X1(D,79),'k--');hold on  %tp50

xlim([Time1 Time2])
ylim([0 15])
set(gca,'xtick',[0:120:1440]+10000)
%set(gca,'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' '])
set(gca,'XGrid','on')
%=========================================================================
%=========================================================================

%Add Model
%=========
SMV2=sbiomodel('SMV2');
%Add species 
%===========
Species={'Ikba','Ikbb','Ikbe','Ikbd','p100','p105','RelA','RelB','p50','p52','R_lB',            ...
         'RelAp50','RelAp52','RelBp50','RelBp52','RelBp100','RelAp100','R_lBp50','R_lBp100',    ...
         'IkbaRelAp50','IkbbRelAp50','IkbeRelAp50','IkbdRelAp50',                               ...
         'IkbaRelAp52','IkbbRelAp52','IkbeRelAp52','IkbdRelAp52',                               ...
         'IkbaRelBp50','IkbbRelBp50','IkbeRelBp50','IkbdRelBp50',                               ...
         'IkbaR_lBp50','IkbbR_lBp50','IkbeR_lBp50','IkbdR_lBp50',                               ...
         'tIkba','tIkbb','tIkbe','tIkbd','tp100','tp105','tRelA','tRelB','tp50','tp52','tR_lB', ...
         'NEMO_ACTIVITY','NIK_ACTIVITY'};
Speciesn=cell(1,34);
for w=1:length(Species)
if w >= 36
continue
end
S=Species{w};S(end+1)='n'; Speciesn{w}=S;
end
I='InitialAmountUnit';U='nanomole';i='InitialAmount';C='ConstantAmount';B='BoundaryCondition';t='true';f='false';v=' ';
for w=1:length(Species)
if w<=35
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
else
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),true,sprintf(i),0,sprintf(I),sprintf(U));']);
end
     if w>=36;continue;end
     eval([sprintf(Speciesn{w}) '=addspecies(SMV2,Speciesn{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
end
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'InitialAmount',1.00);
set (NEMO_ACTIVITY, 'ConstantAmount', false);
set (NEMO_ACTIVITY, 'BoundaryCondition', true);
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'ConstantAmount', false);
set (NIK_ACTIVITY, 'BoundaryCondition', true);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

                       
M='MassAction';n=1;V='ValueUnits'; U0='nanomole/minute';U1='1/minute';U2='(1/nanomole)*(1/minute)';U3='1/nanomole';
for q1=1:11
  if (q1~=8) && (q1~=11)
      %Constitutive transcription
      Rkonst=horzcat('Kct',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst),Parameters.monparams.Const_trnx(q1),V,U0);']);n=n+1; 
  
      %Contitutive mRNA degradation
      Rconst=horzcat('Kmdeg',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.monparams.mRNA_Deg(q1),V,U1);']);n=n+1;
  end
  if (q1~=4) && (q1~=10)
      %Translation Rate
      R3=horzcat(Species{q1+35},' ','->',' ',Species{q1});RO3=horzcat('Trnsln_',Species{q1},num2str(n));KLO3=horzcat('kine',num2str(n));Rkonst3=horzcat('K',num2str(n));
      eval([sprintf(RO3) '=addreaction(SMV2,sprintf(R3));']); 
      eval([sprintf(KLO3) '= addkineticlaw(eval(sprintf(RO3)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst3), Parameters.monparams.Trnsln_rate (q1),V,U1);']);
      set(eval(sprintf(KLO3)),'ParameterVariableName', {sprintf(Rkonst3)});n=n+1;
  end
      %Constitutive Degradation of Free primary Species (Cytoplasm)
      R4=horzcat(Species{q1},' ','->',' ','null');RO4=horzcat('Const_Prot_Deg_',Species{q1},num2str(n));KLO4=horzcat('kine',num2str(n));Rkonst4=horzcat('K',num2str(n));
      eval([sprintf(RO4) '=addreaction(SMV2,sprintf(R4));']); 
      eval([sprintf(KLO4) '= addkineticlaw(eval(sprintf(RO4)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst4),Parameters.monparams.CnstDegFrSpCyt(q1),V,U1);']); 
      set(eval(sprintf(KLO4)),'ParameterVariableName', {sprintf(Rkonst4)});n=n+1;
      %Constitutive Degradation of Free primary Species (Nucleus)
      R5=horzcat(Speciesn{q1},' ','->',' ','null');RO5=horzcat('Const_Prot_Deg_',Speciesn{q1},num2str(n));KLO5=horzcat('kine',num2str(n));Rkonst5=horzcat('K',num2str(n));
      eval([sprintf(RO5) '=addreaction(SMV2,sprintf(R5));']); 
      eval([sprintf(KLO5) '= addkineticlaw(eval(sprintf(RO5)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst5),Parameters.monparams.CnstDegFrSpNuc(q1),V,U1);']); 
      set(eval(sprintf(KLO5)),'ParameterVariableName', {sprintf(Rkonst5)});n=n+1;
      %NEMO Mediated Degradation of Canonical Ikbs
      if q1<=3
      R6=horzcat(Species{q1},' ','+',' ','NEMO_ACTIVITY',' ','->',' ','null');RO6=horzcat('NEMO_Mon_Deg_',Species{q1},num2str(n));KLO6=horzcat('kine',num2str(n));Rkonst6=horzcat('K',num2str(n));
      eval([sprintf(RO6) '=addreaction(SMV2,sprintf(R6));']); 
      eval([sprintf(KLO6) '= addkineticlaw(eval(sprintf(RO6)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst6),Parameters.monparams.NEMODegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO6)),'ParameterVariableName', {sprintf(Rkonst6)});n=n+1;
      end
      %NIK Mediated Degradation of Non-Canonical Ikb
      if q1==4
      R7=horzcat(Species{q1},' ','+',' ','NIK_ACTIVITY',' ','->',' ','null');RO7=horzcat('NIK_Mon_Deg_',Species{q1},num2str(n));KLO7=horzcat('kine',num2str(n));Rkonst7=horzcat('K',num2str(n));
      eval([sprintf(RO7) '=addreaction(SMV2,sprintf(R7));']); 
      eval([sprintf(KLO7) '= addkineticlaw(eval(sprintf(RO7)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst7), Parameters.monparams.NIKDegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO7)),'ParameterVariableName', {sprintf(Rkonst7)});n=n+1;
      end
      %Nuclear Import 
      if q1<=4
      R8=horzcat(Species{q1},' ','->',' ',Speciesn{q1});RO8=horzcat('Nuc_Im_Mon',Species{q1},num2str(n));KLO8=horzcat('kine',num2str(n));Rkonst8=horzcat('K',num2str(n));
      eval([sprintf(RO8) '=addreaction(SMV2,sprintf(R8));']); 
      eval([sprintf(KLO8) '= addkineticlaw(eval(sprintf(RO8)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst8),Parameters.monparams.Nuc_Imp(q1),V,U1);']); 
      set(eval(sprintf(KLO8)),'ParameterVariableName', {sprintf(Rkonst8)});n=n+1;
      end
      %Nuclear Export
      if q1<=4
      R9=horzcat(Speciesn{q1},' ','->',' ',Species{q1});RO9=horzcat('Nuc_Ex_Mon',Species{q1},num2str(n));KLO9=horzcat('kine',num2str(n));Rkonst9=horzcat('K',num2str(n));
      eval([sprintf(RO9) '=addreaction(SMV2,sprintf(R9));']); 
      eval([sprintf(KLO9) '= addkineticlaw(eval(sprintf(RO9)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst9),Parameters.monparams.Nuc_Exp(q1),V,U1);']); 
      set(eval(sprintf(KLO9)),'ParameterVariableName', {sprintf(Rkonst9)});n=n+1;
      end
end
%====================================================1.6e-3 


for q2=1:8
       %Association  of Dimers in Cytoplasm
       R10=horzcat(Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end),' ','->',' ',Species{11+q2});RO10=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO10=horzcat('kine',num2str(n));Rkonst10=horzcat('K',num2str(n));
       eval([sprintf(RO10) '=addreaction(SMV2,sprintf(R10));']); 
       eval([sprintf(KLO10) '= addkineticlaw(eval(sprintf(RO10)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst10),Parameters.dimparams.MonAssCyt(q2),V,U2);']); 
       set(eval(sprintf(KLO10)),'ParameterVariableName', {sprintf(Rkonst10)});n=n+1;
       %Association  of Dimers in Nucleus
       R11=horzcat(Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end),' ','->',' ',Speciesn{11+q2});RO11=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO11=horzcat('kine',num2str(n));Rkonst11=horzcat('K',num2str(n));
       eval([sprintf(RO11) '=addreaction(SMV2,sprintf(R11));']); 
       eval([sprintf(KLO11) '= addkineticlaw(eval(sprintf(RO11)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst11),Parameters.dimparams.MonAssNuc(q2),V,U2);']); 
       set(eval(sprintf(KLO11)),'ParameterVariableName', {sprintf(Rkonst11)});n=n+1;
       %Dissociation  of Dimers in Cytoplasm
       R12=horzcat(Species{11+q2},' ','->',' ',Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end));RO12=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO12=horzcat('kine',num2str(n));Rkonst12=horzcat('K',num2str(n));
       eval([sprintf(RO12) '=addreaction(SMV2,sprintf(R12));']); 
       eval([sprintf(KLO12) '= addkineticlaw(eval(sprintf(RO12)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst12),Parameters.dimparams.DimDissCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO12)),'ParameterVariableName', {sprintf(Rkonst12)});n=n+1;
       %Dissociation  of Dimers in Nucleus
       R13=horzcat(Speciesn{11+q2},' ','->',' ',Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end));RO13=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO13=horzcat('kine',num2str(n));Rkonst13=horzcat('K',num2str(n));
       eval([sprintf(RO13) '=addreaction(SMV2,sprintf(R13));']); 
       eval([sprintf(KLO13) '= addkineticlaw(eval(sprintf(RO13)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst13),Parameters.dimparams.DimDissNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO13)),'ParameterVariableName', {sprintf(Rkonst13)});n=n+1;
       %Nuclear Import
       R14=horzcat(Species{11+q2},' ','->',' ',Speciesn{11+q2});RO14=horzcat('Nuc_Im_Dim',Species{q2+11},num2str(n));KLO14=horzcat('kine',num2str(n));Rkonst14=horzcat('K',num2str(n));
       eval([sprintf(RO14) '=addreaction(SMV2,sprintf(R14));']); 
       eval([sprintf(KLO14) '= addkineticlaw(eval(sprintf(RO14)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst14),Parameters.dimparams.NucImp(q2),V,U1);']); 
       set(eval(sprintf(KLO14)),'ParameterVariableName', {sprintf(Rkonst14)});n=n+1;
       %Nuclear Export
       R15=horzcat(Speciesn{11+q2},' ','->',' ',Species{11+q2});RO15=horzcat('Nuc_Ex_Dim',Species{q2+11},num2str(n));KLO15=horzcat('kine',num2str(n));Rkonst15=horzcat('K',num2str(n));
       eval([sprintf(RO15) '=addreaction(SMV2,sprintf(R15));']); 
       eval([sprintf(KLO15) '= addkineticlaw(eval(sprintf(RO15)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst15),Parameters.dimparams.NucExp(q2),V,U1);']); 
       set(eval(sprintf(KLO15)),'ParameterVariableName', {sprintf(Rkonst15)});n=n+1;
       %Free Dimer Degradation in Cytoplasm
       R16=horzcat(Species{11+q2},' ','->',' ','null');RO16=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO16=horzcat('kine',num2str(n));Rkonst16=horzcat('K',num2str(n));
       eval([sprintf(RO16) '=addreaction(SMV2,sprintf(R16));']); 
       eval([sprintf(KLO16) '= addkineticlaw(eval(sprintf(RO16)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst16),Parameters.dimparams.DimDegCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO16)),'ParameterVariableName', {sprintf(Rkonst16)});n=n+1;
       %Free Dimer Degradation in Nucleus
       R17=horzcat(Speciesn{11+q2},' ','->',' ','null');RO17=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO17=horzcat('kine',num2str(n));Rkonst17=horzcat('K',num2str(n));
       eval([sprintf(RO17) '=addreaction(SMV2,sprintf(R17));']); 
       eval([sprintf(KLO17) '= addkineticlaw(eval(sprintf(RO17)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst17),Parameters.dimparams.DimDegNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO17)),'ParameterVariableName', {sprintf(Rkonst17)});n=n+1;
end
%========================================================
                                                     
for q3=1:16
        %Association of trimers in Cytoplasm
        R18=horzcat(Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end),' ','->',' ',Species{q3+19});RO18=horzcat('Trim_Ass_Cyto',Species{q3+19},num2str(n));KLO18=horzcat('kine',num2str(n));Rkonst18=horzcat('K',num2str(n));
        eval([sprintf(RO18) '=addreaction(SMV2,sprintf(R18));']); 
        eval([sprintf(KLO18) '= addkineticlaw(eval(sprintf(RO18)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst18),Parameters.trimparams.Trim_Ass_Cyt(q3),V,U2);']); 
        set(eval(sprintf(KLO18)),'ParameterVariableName', {sprintf(Rkonst18)});n=n+1;
        %Association of trimers in Nucleus
        R19=horzcat(Speciesn{q3+19}(1:4),'n',' ','+',' ',Speciesn{q3+19}(5:end),' ','->',' ',Speciesn{q3+19});RO19=horzcat('Trim_Ass_Nuc',Species{q3+19},num2str(n));KLO19=horzcat('kine',num2str(n));Rkonst19=horzcat('K',num2str(n));
        eval([sprintf(RO19) '=addreaction(SMV2,sprintf(R19));']); 
        eval([sprintf(KLO19) '= addkineticlaw(eval(sprintf(RO19)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst19),Parameters.trimparams.Trim_Ass_Nuc(q3),V,U2);']); 
        set(eval(sprintf(KLO19)),'ParameterVariableName', {sprintf(Rkonst19)});n=n+1;
        %Dissociation of trimers in Cytoplasm
        R20=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end));RO20=horzcat('Trim_Dis_Cyt',Species{q3+19},num2str(n));KLO20=horzcat('kine',num2str(n));Rkonst20=horzcat('K',num2str(n));
        eval([sprintf(RO20) '=addreaction(SMV2,sprintf(R20));']); 
        eval([sprintf(KLO20) '= addkineticlaw(eval(sprintf(RO20)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst20),Parameters.trimparams.Trim_Dis_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO20)),'ParameterVariableName', {sprintf(Rkonst20)});n=n+1;
        %Dissociation of trimers in Nucleus
        R21=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n',' ','+',' ',Species{q3+19}(5:end),'n');RO21=horzcat('Trim_Dis_Nuc',Species{q3+19},num2str(n));KLO21=horzcat('kine',num2str(n));Rkonst21=horzcat('K',num2str(n));
        eval([sprintf(RO21) '=addreaction(SMV2,sprintf(R21));']); 
        eval([sprintf(KLO21) '= addkineticlaw(eval(sprintf(RO21)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst21),Parameters.trimparams.Trim_Dis_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO21)),'ParameterVariableName', {sprintf(Rkonst21)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Cytoplasm
        R22=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(5:end));RO22=horzcat('Const_Deg_Trim_Ikb_Cyt',Species{q3+19},num2str(n));KLO22=horzcat('kine',num2str(n));Rkonst22=horzcat('K',num2str(n));
        eval([sprintf(RO22) '=addreaction(SMV2,sprintf(R22));']); 
        eval([sprintf(KLO22) '= addkineticlaw(eval(sprintf(RO22)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst22),Parameters.trimparams.Const_Trim_Deg_IkB_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO22)),'ParameterVariableName', {sprintf(Rkonst22)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Nucleus
        R23=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(5:end),'n');RO23=horzcat('Const_Deg_Trim_Ikb_Nuc',Species{q3+19},num2str(n));KLO23=horzcat('kine',num2str(n));Rkonst23=horzcat('K',num2str(n));
        eval([sprintf(RO23) '=addreaction(SMV2,sprintf(R23));']); 
        eval([sprintf(KLO23) '= addkineticlaw(eval(sprintf(RO23)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst23),Parameters.trimparams.Const_Trim_Deg_IkB_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO23)),'ParameterVariableName', {sprintf(Rkonst23)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Cytoplasm
        R24=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4));RO24=horzcat('Const_Deg_Trim_NFkb_Cyt',Species{q3+19},num2str(n));KLO24=horzcat('kine',num2str(n));Rkonst24=horzcat('K',num2str(n));
        eval([sprintf(RO24) '=addreaction(SMV2,sprintf(R24));']); 
        eval([sprintf(KLO24) '= addkineticlaw(eval(sprintf(RO24)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst24),Parameters.trimparams.Const_Trim_Deg_NFkb_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO24)),'ParameterVariableName', {sprintf(Rkonst24)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Nucleus
        R25=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n');RO25=horzcat('Const_Deg_Trim_NFkb_Nuc',Species{q3+19},num2str(n));KLO25=horzcat('kine',num2str(n));Rkonst25=horzcat('K',num2str(n));
        eval([sprintf(RO25) '=addreaction(SMV2,sprintf(R25));']); 
        eval([sprintf(KLO25) '= addkineticlaw(eval(sprintf(RO25)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst25),Parameters.trimparams.Const_Trim_Deg_NFkb_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO25)),'ParameterVariableName', {sprintf(Rkonst25)});n=n+1;
        %NEMO Mediated Degradation of Ikb in Trimer in Cytoplasm
        R26=horzcat(Species{q3+19},' ','+',' ','NEMO_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO26=horzcat('NEMO_DEg_inTrim',Species{q3+19},num2str(n));KLO26=horzcat('kine',num2str(n));Rkonst26=horzcat('K',num2str(n));
        eval([sprintf(RO26) '=addreaction(SMV2,sprintf(R26));']); 
        eval([sprintf(KLO26) '= addkineticlaw(eval(sprintf(RO26)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst26),Parameters.trimparams.NEMO_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO26)),'ParameterVariableName', {sprintf(Rkonst26)});n=n+1;
        %NIK Mediated Degradation of Ikb in Trimer in Cytoplasm
        R27=horzcat(Species{q3+19},' ','+',' ','NIK_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO27=horzcat('NIK_DEg_inTrim',Species{q3+19},num2str(n));KLO27=horzcat('kine',num2str(n));Rkonst27=horzcat('K',num2str(n));
        eval([sprintf(RO27) '=addreaction(SMV2,sprintf(R27));']); 
        eval([sprintf(KLO27) '= addkineticlaw(eval(sprintf(RO27)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst27),Parameters.trimparams.NIK_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO27)),'ParameterVariableName', {sprintf(Rkonst27)});n=n+1;
        %Nuclear Import of Trimers
        R28=horzcat(Species{q3+18},' ','->',Speciesn{q3+18});RO28=horzcat('Nuc_Imp_Trim',Species{q3+18},num2str(n));KLO28=horzcat('kine',num2str(n));Rkonst28=horzcat('K',num2str(n));
        eval([sprintf(RO28) '=addreaction(SMV2,sprintf(R28));']); 
        eval([sprintf(KLO28) '= addkineticlaw(eval(sprintf(RO28)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst28),Parameters.trimparams.Nuc_Imp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO28)),'ParameterVariableName', {sprintf(Rkonst28)});n=n+1;
        %Nuclear Export of Trimers
        R29=horzcat(Speciesn{q3+19},' ','->',Species{q3+19});RO29=horzcat('Nuc_Exp_Trip',Species{q3+19},num2str(n));KLO29=horzcat('kine',num2str(n));Rkonst29=horzcat('K',num2str(n));
        eval([sprintf(RO29) '=addreaction(SMV2,sprintf(R29));']); 
        eval([sprintf(KLO29) '= addkineticlaw(eval(sprintf(RO29)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst29),Parameters.trimparams.Nuc_Exp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO29)),'ParameterVariableName', {sprintf(Rkonst29)});n=n+1;
end

   
        R30=horzcat('2',' ',Species{5},' ','->',Species{4});RO30=horzcat('p100_Ass',Species{5},num2str(n));KLO30=horzcat('kine',num2str(n));Rkonst30=horzcat('K',num2str(n));
        eval([sprintf(RO30) '=addreaction(SMV2,sprintf(R30));']); 
        eval([sprintf(KLO30) '= addkineticlaw(eval(sprintf(RO30)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst30),Parameters.p100_Ass,V,U2);']); 
        set(eval(sprintf(KLO30)),'ParameterVariableName', {sprintf(Rkonst30)});n=n+1;

        R31=horzcat(Species{4},' ','->',' ','2',' ',Species{5});RO31=horzcat('p100_Ass',Species{4},num2str(n));KLO31=horzcat('kine',num2str(n));Rkonst31=horzcat('K',num2str(n));
        eval([sprintf(RO31) '=addreaction(SMV2,sprintf(R31));']); 
        eval([sprintf(KLO31) '= addkineticlaw(eval(sprintf(RO31)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst31),Parameters.Ikbd_Diss,V,U1);']); 
        set(eval(sprintf(KLO31)),'ParameterVariableName', {sprintf(Rkonst31)});n=n+1;

        R32=horzcat(Species{5},' ','+',' ',Species{end},' ','->',' ',Species{10});RO32=horzcat('p52_processing',Species{5},num2str(n));KLO32=horzcat('kine',num2str(n));Rkonst32=horzcat('K',num2str(n));
        eval([sprintf(RO32) '=addreaction(SMV2,sprintf(R32));']); 
        eval([sprintf(KLO32) '= addkineticlaw(eval(sprintf(RO32)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst32),Parameters.NIK_mediated_p52,V,U2);']); 
        set(eval(sprintf(KLO32)),'ParameterVariableName', {sprintf(Rkonst32)});n=n+1;
%===========================================================
%Chequer Boxes of Inducibility 
%=============================
%                   
%Add Rules
%==========

RfntIkba = 'tIkba = ( Kct1 * ( ( 1 + ( W11 * ( RelAp50n / Kd11 ) ^ 3 ) + ( W12 * ( RelAp52n / Kd12 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd11 ) ^ 3 + ( RelAp52n / Kd12 ) ^ 3 ) ) ) - Kmdeg1 * tIkba';
              
RfntIkbb = 'tIkbb = ( Kct2 * ( ( 1 + ( W21 * ( RelAp50n / Kd21 ) ^ 3 ) + ( W22 * ( RelAp52n / Kd22 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd21 ) ^ 3 + ( RelAp52n / Kd22 ) ^ 3 ) ) ^ Kdel1 ) - Kmdeg3 * tIkbb';
                 
RfntIkbe = 'tIkbe = ( Kct3 * ( ( 1 + ( W31 * ( RelAp50n / Kd31 ) ^ 3 ) + ( W32 * ( RelAp52n / Kd32 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd31 ) ^ 3 + ( RelAp52n / Kd32 ) ^ 3 ) ) ^ Kdel2 ) - Kmdeg3 * tIkbe';
                
Rfntp100 = 'tp100 = ( Kct5 * ( ( 1 + ( W41 * ( RelAp50n / Kd41 ) ^ 3 ) + ( W42 * ( RelAp52n / Kd42 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd41 ) ^ 3 + ( RelAp52n / Kd42 ) ^ 3 ) ) ^ Kdel3 ) - Kmdeg5 * tp100';
                 
RfntRelA = 'tRelA = Kct7  - Kmdeg7 * tRelA';

RfntRelB = 'tRelB = ( Kct8 * ( ( 1 + ( W51 * ( RelAp50n * Kdel5 / Kd51 ) ^ 3 ) + ( W52 * ( RelAp52n * Kdel5 / Kd52 ) ^ 3 ) + ( W53 * ( RelBp50n / Kd53 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd51 ) ^ 3 + ( RelAp52n * Kdel5 / Kd52 ) ^ 3  + ( RelBp50n / Kd53 ) ^ 3 ) ) ) - Kmdeg8 * tRelB';

%RfntRelB = 'tRelB = ( Kct8 * ( ( 1 + ( W51 * ( RelAp50n / Kd51 ) ^ 3 ) + ( W52 * ( RelAp52n / Kd52 ) ^ 3 ) + ( W53 * ( RelBp50n / Kd53 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd51 ) ^ 3 + ( RelAp52n / Kd52 ) ^ 3  + ( RelBp50n / Kd53 ) ^ 3 ) ) ^ Kdel4 ) - Kmdeg8 * tRelB';

RfntR_lB = 'tR_lB = ( Kct11 * ( ( 1 + ( W61 * ( RelAp50n * Kdel5 / Kd61 ) ^ 3 ) + ( W62 * ( RelAp52n * Kdel5 / Kd62 ) ^ 3 ) + ( W63 * ( RelBp50n / Kd63 ) ^ 3 ) + ( W64 * ( R_lBp50n / Kd64 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd61 ) ^ 3 + ( RelAp52n * Kdel5 / Kd62 ) ^ 3  + ( RelBp50n / Kd63 ) ^ 3 + ( R_lBp50n / Kd64 ) ^ 3 ) ) ^ Kdel4 ) - Kmdeg11 * tR_lB';
                   
Rfntp50  = 'tp50  = Kct9 - Kmdeg9 * tp50';
for krow= 1:size(Parameters.Kd,1)
  for kcolumn= 1:size(Parameters.Kd,2)
      Rconst=horzcat('Kd',num2str(krow),num2str(kcolumn)); Rconst1=horzcat('W',num2str(krow),num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.Kd(krow,kcolumn),V,U);']);n=n+1;
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst1),Parameters.W(krow,kcolumn));']);n=n+1; 
  end
end
%============================================
  RlO1 = addrule(SMV2,sprintf(RfntIkba),'rate');
  RlO2 = addrule(SMV2,sprintf(RfntIkbb),'rate');
  RlO3 = addrule(SMV2,sprintf(RfntIkbe),'rate');
  RlO4 = addrule(SMV2,sprintf(Rfntp100),'rate');
  RlO5 = addrule(SMV2,sprintf(RfntRelA),'rate');
  RlO6 = addrule(SMV2,sprintf(RfntRelB),'rate');
  RlO7 = addrule(SMV2,sprintf(RfntR_lB),'rate');
  RlO8 = addrule(SMV2,sprintf(Rfntp50),'rate');
%============================================
Ca='ConstantValue';%VA='Value'

%DELAY
%=====
%kdels are 1 during euilibration
for kcolumn= 1:5
      Rconst=horzcat('Kdel',num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),1,Ca,false);']);
      n=n+1; 
end
%kdels are 0 at stimulation 0 min.
for P=1:4
 E=horzcat('Eee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(0));
                                         T=horzcat('time >= ',num2str(10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%kdels are 1 at stimulation phase once delay is over.
DELAY=[37 37 45 0]; %beta epsilon p100 RelB 
for P=1:4
 E=horzcat('Eeee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(1));
                                         T=horzcat('time >= ',num2str(DELAY(P)+10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%RelB star issue
%===============

      RconsT=horzcat('Kct',num2str(8));RconsTd=horzcat('Kmdeg',num2str(8));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsT),Parameters.monparams.Const_trnx(8),Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsTd),Parameters.monparams.mRNA_Deg(8),Ca,false);']); n=n+1; 
      EEE1=addevent(SMV2,'time >= 10000','Kct8 = 0');%EEE2=addevent(SMV2,'time >= 10000','Kmdeg8 = 0');

      RconsTT=horzcat('Kct',num2str(11));RconsTTd=horzcat('Kmdeg',num2str(11));
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTT),0,Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTTd),0,Ca,false);']); n=n+1; 
      %EEE3=addevent(SMV2,'time >= 150000','Kct11 =  3.6e-5');EEE4=addevent(SMV2,'time >= 150000','Kmdeg11 = 4.6e-3');
      E=horzcat('EEE',num2str(3));N=horzcat('Kct',num2str(11),' ','=',' ',num2str(Parameters.monparams.Const_trnx(11))); T=horzcat('time >= ',num2str(10000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);       
      E=horzcat('EEE',num2str(4));N=horzcat('Kmdeg',num2str(11),' ','=',' ',num2str(Parameters.monparams.mRNA_Deg(11))); %T=horzcat('time >= ',num2str(150000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);                                                                                                

%Effective RelA heterodimer concentration in nucleus involved with RelB synthesis
evalues =  [1 .01];etime = [0 60];
 for Q=0:1:60
 Kdel5 = interp1(etime, evalues, Q,'pchip');E=horzcat('Eeff',num2str(Q));N=horzcat('Kdel5',' ','=',' ',num2str(Kdel5));
                                         T=horzcat('time >= ',num2str(10030+Q));%RelA effective concentration start going down after 30 mins
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end     
      
%SIGNAL
%======
%TNFc
%cvalues = [1 60 100 65 50 36 20 10 3];ctime = [0 5 10 15 20 25 30 360 480];
cvalues =  [1 60 100 65 50 36 20 18 15];ctime = [0 5 10 15 20 25 30 360 480];
%======================================
%Add Parameters of Inducibility
%==============================
for Q=0:1:480
 NEMO = interp1(ctime, cvalues, Q,'pchip');E=horzcat('E',num2str(Q));N=horzcat('NEMO_ACTIVITY',' ','=',' ',num2str(NEMO));
                                         T=horzcat('time >= ',num2str(10000+Q));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end

%SIMULATE MODEL
%==============
CnfgstO = addconfigset(SMV2, 'De3');
set(CnfgstO.CompileOptions, 'UnitConversion', false);
set(CnfgstO.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO,'TimeUnits',' ');
set(CnfgstO,'SolverType','sundials');
set(CnfgstO,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO)
[t1,x1,names1]= sbiosimulate(SMV2,CnfgstO);
%ADD VARIANT 1
%==============
estvarObj1 = addvariant (SMV2, 'nfkb2kn');
addcontent(estvarObj1, {'parameter','K34', 'Value',0});
CnfgstO1 = addconfigset(SMV2, 'nfkb2kn');
set(CnfgstO1.CompileOptions, 'UnitConversion', false);
set(CnfgstO1.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO1,'TimeUnits',' ');
set(CnfgstO1,'SolverType','sundials');
set(CnfgstO1,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO1,estvarObj1)
[t2,x2,names2]= sbiosimulate(SMV2,CnfgstO1,estvarObj1);
%ADD VARIANT 2
%==============
estvarObj2 = addvariant (SMV2, 'RelAnfkb2kn');
addcontent(estvarObj2, {'parameter','K34', 'Value',0});
addcontent(estvarObj2, {'parameter','K44', 'Value',0});
%addcontent(estvarObj2, {'parameter','Kct1', 'Value',0.5e-3});
CnfgstO2 = addconfigset(SMV2, 'RelAnfkb2kn');
set(CnfgstO2.CompileOptions, 'UnitConversion', false);
set(CnfgstO2.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO2,'TimeUnits',' ');
set(CnfgstO2,'SolverType','sundials');
set(CnfgstO2,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO2,estvarObj2)
[t3,x3,names3]= sbiosimulate(SMV2,CnfgstO2,estvarObj2);

%COMMIT MODEL TO PULSE
%=====================

 pvalues = [1 60 100 65 50 36 21 16 10 1 1];ptime = [0 5 10 15 20 25 30 45 59 60 480];
for Q=0:1:480
 NEMO = interp1(ptime, pvalues, Q,'pchip');E=horzcat('E',num2str(Q));N=horzcat('NEMO_ACTIVITY',' ','=',' ',num2str(NEMO));
                                           T=horzcat('time >= ',num2str(10000+Q));
 %eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
 set(eval(sprintf(E)),'EventFcns', {sprintf(N)});
end

sbioaccelerate(SMV2,CnfgstO)
[t7,x7,names4]= sbiosimulate(SMV2,CnfgstO);

%ADD VARIANT 3
%==============
estvarObj3 = addvariant (SMV2, 'nfkb2knPulse');
addcontent(estvarObj3, {'parameter','K34', 'Value',0});


CnfgstO3 = addconfigset(SMV2, 'nfkb2knPulse');
set(CnfgstO3.CompileOptions, 'UnitConversion', false);
set(CnfgstO3.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO3,'TimeUnits',' ');
set(CnfgstO3,'SolverType','sundials');
set(CnfgstO3,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO3,estvarObj3)
[t8,x8,names5]= sbiosimulate(SMV2,CnfgstO3,estvarObj3);

%ADD VARIANT 4
%==============
estvarObj2 = addvariant (SMV2, 'RelAnfkb2knPulse');
addcontent(estvarObj2, {'parameter','K34', 'Value',0});
addcontent(estvarObj2, {'parameter','K44', 'Value',0});
%addcontent(estvarObj2, {'parameter','Kct1', 'Value',0.5e-3});
CnfgstO2 = addconfigset(SMV2, 'RelAnfkb2knPulse');
set(CnfgstO2.CompileOptions, 'UnitConversion', false);
set(CnfgstO2.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO2,'TimeUnits',' ');
set(CnfgstO2,'SolverType','sundials');
set(CnfgstO2,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO2,estvarObj2)
[t9,x9,names3]= sbiosimulate(SMV2,CnfgstO2,estvarObj2);


%Add Model
%=========
SMV2=sbiomodel('SMV2');
%Add species 
%===========
Species={'Ikba','Ikbb','Ikbe','Ikbd','p100','p105','RelA','RelB','p50','p52','R_lB',            ...
         'RelAp50','RelAp52','RelBp50','RelBp52','RelBp100','RelAp100','R_lBp50','R_lBp100',    ...
         'IkbaRelAp50','IkbbRelAp50','IkbeRelAp50','IkbdRelAp50',                               ...
         'IkbaRelAp52','IkbbRelAp52','IkbeRelAp52','IkbdRelAp52',                               ...
         'IkbaRelBp50','IkbbRelBp50','IkbeRelBp50','IkbdRelBp50',                               ...
         'IkbaR_lBp50','IkbbR_lBp50','IkbeR_lBp50','IkbdR_lBp50',                               ...
         'tIkba','tIkbb','tIkbe','tIkbd','tp100','tp105','tRelA','tRelB','tp50','tp52','tR_lB', ...
         'NEMO_ACTIVITY','NIK_ACTIVITY'};
Speciesn=cell(1,34);
for w=1:length(Species)
if w >= 36
continue
end
S=Species{w};S(end+1)='n'; Speciesn{w}=S;
end
I='InitialAmountUnit';U='nanomole';i='InitialAmount';C='ConstantAmount';B='BoundaryCondition';t='true';f='false';v=' ';
for w=1:length(Species)
if w<=35
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
else
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),true,sprintf(i),0,sprintf(I),sprintf(U));']);
end
     if w>=36;continue;end
     eval([sprintf(Speciesn{w}) '=addspecies(SMV2,Speciesn{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
end
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

set (NEMO_ACTIVITY, 'ConstantAmount', false);
set (NEMO_ACTIVITY, 'BoundaryCondition', true);
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'ConstantAmount', false);
set (NIK_ACTIVITY, 'BoundaryCondition', true);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

Parameters.monparams.Trnsln_rate   =  [12           12         12         0          0          0         12         12          12          0        12        ];
    
for q1=1:11
  if (q1~=8) && (q1~=11)
      %Constitutive transcription
      Rkonst=horzcat('Kct',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst),Parameters.monparams.Const_trnx(q1),V,U0);']);n=n+1; 
  
      %Contitutive mRNA degradation
      Rconst=horzcat('Kmdeg',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.monparams.mRNA_Deg(q1),V,U1);']);n=n+1;
  end
  if (q1~=4) && (q1~=10)
      %Translation Rate
      R3=horzcat(Species{q1+35},' ','->',' ',Species{q1});RO3=horzcat('Trnsln_',Species{q1},num2str(n));KLO3=horzcat('kine',num2str(n));Rkonst3=horzcat('K',num2str(n));
      eval([sprintf(RO3) '=addreaction(SMV2,sprintf(R3));']); 
      eval([sprintf(KLO3) '= addkineticlaw(eval(sprintf(RO3)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst3), Parameters.monparams.Trnsln_rate (q1),V,U1);']);
      set(eval(sprintf(KLO3)),'ParameterVariableName', {sprintf(Rkonst3)});n=n+1;
  end
      %Constitutive Degradation of Free primary Species (Cytoplasm)
      R4=horzcat(Species{q1},' ','->',' ','null');RO4=horzcat('Const_Prot_Deg_',Species{q1},num2str(n));KLO4=horzcat('kine',num2str(n));Rkonst4=horzcat('K',num2str(n));
      eval([sprintf(RO4) '=addreaction(SMV2,sprintf(R4));']); 
      eval([sprintf(KLO4) '= addkineticlaw(eval(sprintf(RO4)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst4),Parameters.monparams.CnstDegFrSpCyt(q1),V,U1);']); 
      set(eval(sprintf(KLO4)),'ParameterVariableName', {sprintf(Rkonst4)});n=n+1;
      %Constitutive Degradation of Free primary Species (Nucleus)
      R5=horzcat(Speciesn{q1},' ','->',' ','null');RO5=horzcat('Const_Prot_Deg_',Speciesn{q1},num2str(n));KLO5=horzcat('kine',num2str(n));Rkonst5=horzcat('K',num2str(n));
      eval([sprintf(RO5) '=addreaction(SMV2,sprintf(R5));']); 
      eval([sprintf(KLO5) '= addkineticlaw(eval(sprintf(RO5)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst5),Parameters.monparams.CnstDegFrSpNuc(q1),V,U1);']); 
      set(eval(sprintf(KLO5)),'ParameterVariableName', {sprintf(Rkonst5)});n=n+1;
      %NEMO Mediated Degradation of Canonical Ikbs
      if q1<=3
      R6=horzcat(Species{q1},' ','+',' ','NEMO_ACTIVITY',' ','->',' ','null');RO6=horzcat('NEMO_Mon_Deg_',Species{q1},num2str(n));KLO6=horzcat('kine',num2str(n));Rkonst6=horzcat('K',num2str(n));
      eval([sprintf(RO6) '=addreaction(SMV2,sprintf(R6));']); 
      eval([sprintf(KLO6) '= addkineticlaw(eval(sprintf(RO6)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst6),Parameters.monparams.NEMODegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO6)),'ParameterVariableName', {sprintf(Rkonst6)});n=n+1;
      end
      %NIK Mediated Degradation of Non-Canonical Ikb
      if q1==4
      R7=horzcat(Species{q1},' ','+',' ','NIK_ACTIVITY',' ','->',' ','null');RO7=horzcat('NIK_Mon_Deg_',Species{q1},num2str(n));KLO7=horzcat('kine',num2str(n));Rkonst7=horzcat('K',num2str(n));
      eval([sprintf(RO7) '=addreaction(SMV2,sprintf(R7));']); 
      eval([sprintf(KLO7) '= addkineticlaw(eval(sprintf(RO7)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst7), Parameters.monparams.NIKDegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO7)),'ParameterVariableName', {sprintf(Rkonst7)});n=n+1;
      end
      %Nuclear Import 
      if q1<=4
      R8=horzcat(Species{q1},' ','->',' ',Speciesn{q1});RO8=horzcat('Nuc_Im_Mon',Species{q1},num2str(n));KLO8=horzcat('kine',num2str(n));Rkonst8=horzcat('K',num2str(n));
      eval([sprintf(RO8) '=addreaction(SMV2,sprintf(R8));']); 
      eval([sprintf(KLO8) '= addkineticlaw(eval(sprintf(RO8)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst8),Parameters.monparams.Nuc_Imp(q1),V,U1);']); 
      set(eval(sprintf(KLO8)),'ParameterVariableName', {sprintf(Rkonst8)});n=n+1;
      end
      %Nuclear Export
      if q1<=4
      R9=horzcat(Speciesn{q1},' ','->',' ',Species{q1});RO9=horzcat('Nuc_Ex_Mon',Species{q1},num2str(n));KLO9=horzcat('kine',num2str(n));Rkonst9=horzcat('K',num2str(n));
      eval([sprintf(RO9) '=addreaction(SMV2,sprintf(R9));']); 
      eval([sprintf(KLO9) '= addkineticlaw(eval(sprintf(RO9)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst9),Parameters.monparams.Nuc_Exp(q1),V,U1);']); 
      set(eval(sprintf(KLO9)),'ParameterVariableName', {sprintf(Rkonst9)});n=n+1;
      end
end
%====================================================1.6e-3 


for q2=1:8
       %Association  of Dimers in Cytoplasm
       R10=horzcat(Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end),' ','->',' ',Species{11+q2});RO10=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO10=horzcat('kine',num2str(n));Rkonst10=horzcat('K',num2str(n));
       eval([sprintf(RO10) '=addreaction(SMV2,sprintf(R10));']); 
       eval([sprintf(KLO10) '= addkineticlaw(eval(sprintf(RO10)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst10),Parameters.dimparams.MonAssCyt(q2),V,U2);']); 
       set(eval(sprintf(KLO10)),'ParameterVariableName', {sprintf(Rkonst10)});n=n+1;
       %Association  of Dimers in Nucleus
       R11=horzcat(Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end),' ','->',' ',Speciesn{11+q2});RO11=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO11=horzcat('kine',num2str(n));Rkonst11=horzcat('K',num2str(n));
       eval([sprintf(RO11) '=addreaction(SMV2,sprintf(R11));']); 
       eval([sprintf(KLO11) '= addkineticlaw(eval(sprintf(RO11)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst11),Parameters.dimparams.MonAssNuc(q2),V,U2);']); 
       set(eval(sprintf(KLO11)),'ParameterVariableName', {sprintf(Rkonst11)});n=n+1;
       %Dissociation  of Dimers in Cytoplasm
       R12=horzcat(Species{11+q2},' ','->',' ',Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end));RO12=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO12=horzcat('kine',num2str(n));Rkonst12=horzcat('K',num2str(n));
       eval([sprintf(RO12) '=addreaction(SMV2,sprintf(R12));']); 
       eval([sprintf(KLO12) '= addkineticlaw(eval(sprintf(RO12)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst12),Parameters.dimparams.DimDissCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO12)),'ParameterVariableName', {sprintf(Rkonst12)});n=n+1;
       %Dissociation  of Dimers in Nucleus
       R13=horzcat(Speciesn{11+q2},' ','->',' ',Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end));RO13=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO13=horzcat('kine',num2str(n));Rkonst13=horzcat('K',num2str(n));
       eval([sprintf(RO13) '=addreaction(SMV2,sprintf(R13));']); 
       eval([sprintf(KLO13) '= addkineticlaw(eval(sprintf(RO13)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst13),Parameters.dimparams.DimDissNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO13)),'ParameterVariableName', {sprintf(Rkonst13)});n=n+1;
       %Nuclear Import
       R14=horzcat(Species{11+q2},' ','->',' ',Speciesn{11+q2});RO14=horzcat('Nuc_Im_Dim',Species{q2+11},num2str(n));KLO14=horzcat('kine',num2str(n));Rkonst14=horzcat('K',num2str(n));
       eval([sprintf(RO14) '=addreaction(SMV2,sprintf(R14));']); 
       eval([sprintf(KLO14) '= addkineticlaw(eval(sprintf(RO14)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst14),Parameters.dimparams.NucImp(q2),V,U1);']); 
       set(eval(sprintf(KLO14)),'ParameterVariableName', {sprintf(Rkonst14)});n=n+1;
       %Nuclear Export
       R15=horzcat(Speciesn{11+q2},' ','->',' ',Species{11+q2});RO15=horzcat('Nuc_Ex_Dim',Species{q2+11},num2str(n));KLO15=horzcat('kine',num2str(n));Rkonst15=horzcat('K',num2str(n));
       eval([sprintf(RO15) '=addreaction(SMV2,sprintf(R15));']); 
       eval([sprintf(KLO15) '= addkineticlaw(eval(sprintf(RO15)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst15),Parameters.dimparams.NucExp(q2),V,U1);']); 
       set(eval(sprintf(KLO15)),'ParameterVariableName', {sprintf(Rkonst15)});n=n+1;
       %Free Dimer Degradation in Cytoplasm
       R16=horzcat(Species{11+q2},' ','->',' ','null');RO16=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO16=horzcat('kine',num2str(n));Rkonst16=horzcat('K',num2str(n));
       eval([sprintf(RO16) '=addreaction(SMV2,sprintf(R16));']); 
       eval([sprintf(KLO16) '= addkineticlaw(eval(sprintf(RO16)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst16),Parameters.dimparams.DimDegCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO16)),'ParameterVariableName', {sprintf(Rkonst16)});n=n+1;
       %Free Dimer Degradation in Nucleus
       R17=horzcat(Speciesn{11+q2},' ','->',' ','null');RO17=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO17=horzcat('kine',num2str(n));Rkonst17=horzcat('K',num2str(n));
       eval([sprintf(RO17) '=addreaction(SMV2,sprintf(R17));']); 
       eval([sprintf(KLO17) '= addkineticlaw(eval(sprintf(RO17)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst17),Parameters.dimparams.DimDegNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO17)),'ParameterVariableName', {sprintf(Rkonst17)});n=n+1;
end
%========================================================
                                                     
for q3=1:16
        %Association of trimers in Cytoplasm
        R18=horzcat(Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end),' ','->',' ',Species{q3+19});RO18=horzcat('Trim_Ass_Cyto',Species{q3+19},num2str(n));KLO18=horzcat('kine',num2str(n));Rkonst18=horzcat('K',num2str(n));
        eval([sprintf(RO18) '=addreaction(SMV2,sprintf(R18));']); 
        eval([sprintf(KLO18) '= addkineticlaw(eval(sprintf(RO18)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst18),Parameters.trimparams.Trim_Ass_Cyt(q3),V,U2);']); 
        set(eval(sprintf(KLO18)),'ParameterVariableName', {sprintf(Rkonst18)});n=n+1;
        %Association of trimers in Nucleus
        R19=horzcat(Speciesn{q3+19}(1:4),'n',' ','+',' ',Speciesn{q3+19}(5:end),' ','->',' ',Speciesn{q3+19});RO19=horzcat('Trim_Ass_Nuc',Species{q3+19},num2str(n));KLO19=horzcat('kine',num2str(n));Rkonst19=horzcat('K',num2str(n));
        eval([sprintf(RO19) '=addreaction(SMV2,sprintf(R19));']); 
        eval([sprintf(KLO19) '= addkineticlaw(eval(sprintf(RO19)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst19),Parameters.trimparams.Trim_Ass_Nuc(q3),V,U2);']); 
        set(eval(sprintf(KLO19)),'ParameterVariableName', {sprintf(Rkonst19)});n=n+1;
        %Dissociation of trimers in Cytoplasm
        R20=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end));RO20=horzcat('Trim_Dis_Cyt',Species{q3+19},num2str(n));KLO20=horzcat('kine',num2str(n));Rkonst20=horzcat('K',num2str(n));
        eval([sprintf(RO20) '=addreaction(SMV2,sprintf(R20));']); 
        eval([sprintf(KLO20) '= addkineticlaw(eval(sprintf(RO20)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst20),Parameters.trimparams.Trim_Dis_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO20)),'ParameterVariableName', {sprintf(Rkonst20)});n=n+1;
        %Dissociation of trimers in Nucleus
        R21=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n',' ','+',' ',Species{q3+19}(5:end),'n');RO21=horzcat('Trim_Dis_Nuc',Species{q3+19},num2str(n));KLO21=horzcat('kine',num2str(n));Rkonst21=horzcat('K',num2str(n));
        eval([sprintf(RO21) '=addreaction(SMV2,sprintf(R21));']); 
        eval([sprintf(KLO21) '= addkineticlaw(eval(sprintf(RO21)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst21),Parameters.trimparams.Trim_Dis_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO21)),'ParameterVariableName', {sprintf(Rkonst21)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Cytoplasm
        R22=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(5:end));RO22=horzcat('Const_Deg_Trim_Ikb_Cyt',Species{q3+19},num2str(n));KLO22=horzcat('kine',num2str(n));Rkonst22=horzcat('K',num2str(n));
        eval([sprintf(RO22) '=addreaction(SMV2,sprintf(R22));']); 
        eval([sprintf(KLO22) '= addkineticlaw(eval(sprintf(RO22)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst22),Parameters.trimparams.Const_Trim_Deg_IkB_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO22)),'ParameterVariableName', {sprintf(Rkonst22)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Nucleus
        R23=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(5:end),'n');RO23=horzcat('Const_Deg_Trim_Ikb_Nuc',Species{q3+19},num2str(n));KLO23=horzcat('kine',num2str(n));Rkonst23=horzcat('K',num2str(n));
        eval([sprintf(RO23) '=addreaction(SMV2,sprintf(R23));']); 
        eval([sprintf(KLO23) '= addkineticlaw(eval(sprintf(RO23)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst23),Parameters.trimparams.Const_Trim_Deg_IkB_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO23)),'ParameterVariableName', {sprintf(Rkonst23)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Cytoplasm
        R24=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4));RO24=horzcat('Const_Deg_Trim_NFkb_Cyt',Species{q3+19},num2str(n));KLO24=horzcat('kine',num2str(n));Rkonst24=horzcat('K',num2str(n));
        eval([sprintf(RO24) '=addreaction(SMV2,sprintf(R24));']); 
        eval([sprintf(KLO24) '= addkineticlaw(eval(sprintf(RO24)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst24),Parameters.trimparams.Const_Trim_Deg_NFkb_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO24)),'ParameterVariableName', {sprintf(Rkonst24)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Nucleus
        R25=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n');RO25=horzcat('Const_Deg_Trim_NFkb_Nuc',Species{q3+19},num2str(n));KLO25=horzcat('kine',num2str(n));Rkonst25=horzcat('K',num2str(n));
        eval([sprintf(RO25) '=addreaction(SMV2,sprintf(R25));']); 
        eval([sprintf(KLO25) '= addkineticlaw(eval(sprintf(RO25)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst25),Parameters.trimparams.Const_Trim_Deg_NFkb_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO25)),'ParameterVariableName', {sprintf(Rkonst25)});n=n+1;
        %NEMO Mediated Degradation of Ikb in Trimer in Cytoplasm
        R26=horzcat(Species{q3+19},' ','+',' ','NEMO_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO26=horzcat('NEMO_DEg_inTrim',Species{q3+19},num2str(n));KLO26=horzcat('kine',num2str(n));Rkonst26=horzcat('K',num2str(n));
        eval([sprintf(RO26) '=addreaction(SMV2,sprintf(R26));']); 
        eval([sprintf(KLO26) '= addkineticlaw(eval(sprintf(RO26)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst26),Parameters.trimparams.NEMO_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO26)),'ParameterVariableName', {sprintf(Rkonst26)});n=n+1;
        %NIK Mediated Degradation of Ikb in Trimer in Cytoplasm
        R27=horzcat(Species{q3+19},' ','+',' ','NIK_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO27=horzcat('NIK_DEg_inTrim',Species{q3+19},num2str(n));KLO27=horzcat('kine',num2str(n));Rkonst27=horzcat('K',num2str(n));
        eval([sprintf(RO27) '=addreaction(SMV2,sprintf(R27));']); 
        eval([sprintf(KLO27) '= addkineticlaw(eval(sprintf(RO27)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst27),Parameters.trimparams.NIK_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO27)),'ParameterVariableName', {sprintf(Rkonst27)});n=n+1;
        %Nuclear Import of Trimers
        R28=horzcat(Species{q3+18},' ','->',Speciesn{q3+18});RO28=horzcat('Nuc_Imp_Trim',Species{q3+18},num2str(n));KLO28=horzcat('kine',num2str(n));Rkonst28=horzcat('K',num2str(n));
        eval([sprintf(RO28) '=addreaction(SMV2,sprintf(R28));']); 
        eval([sprintf(KLO28) '= addkineticlaw(eval(sprintf(RO28)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst28),Parameters.trimparams.Nuc_Imp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO28)),'ParameterVariableName', {sprintf(Rkonst28)});n=n+1;
        %Nuclear Export of Trimers
        R29=horzcat(Speciesn{q3+19},' ','->',Species{q3+19});RO29=horzcat('Nuc_Exp_Trip',Species{q3+19},num2str(n));KLO29=horzcat('kine',num2str(n));Rkonst29=horzcat('K',num2str(n));
        eval([sprintf(RO29) '=addreaction(SMV2,sprintf(R29));']); 
        eval([sprintf(KLO29) '= addkineticlaw(eval(sprintf(RO29)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst29),Parameters.trimparams.Nuc_Exp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO29)),'ParameterVariableName', {sprintf(Rkonst29)});n=n+1;
end
%=========================================================
%Special Reactions
   Parameters.p100_Ass            = 9e-2;     % p100 association rate to form IkBd in cytoplasm %%% 1.2E-3
   %                           K253 
   Parameters.Ikbd_Diss           = 1.2e-5;   % IkBd dissocation rate to form p100 in cytoplasm %%% 1.2E-2
   %                           K254
   Parameters.NIK_mediated_p52    = 3.5e-3;   % NIK mediated generation of P52 through processing
   %                           K255
   
        R30=horzcat('2',' ',Species{5},' ','->',Species{4});RO30=horzcat('p100_Ass',Species{5},num2str(n));KLO30=horzcat('kine',num2str(n));Rkonst30=horzcat('K',num2str(n));
        eval([sprintf(RO30) '=addreaction(SMV2,sprintf(R30));']); 
        eval([sprintf(KLO30) '= addkineticlaw(eval(sprintf(RO30)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst30),Parameters.p100_Ass,V,U2);']); 
        set(eval(sprintf(KLO30)),'ParameterVariableName', {sprintf(Rkonst30)});n=n+1;

        R31=horzcat(Species{4},' ','->',' ','2',' ',Species{5});RO31=horzcat('p100_Ass',Species{4},num2str(n));KLO31=horzcat('kine',num2str(n));Rkonst31=horzcat('K',num2str(n));
        eval([sprintf(RO31) '=addreaction(SMV2,sprintf(R31));']); 
        eval([sprintf(KLO31) '= addkineticlaw(eval(sprintf(RO31)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst31),Parameters.Ikbd_Diss,V,U1);']); 
        set(eval(sprintf(KLO31)),'ParameterVariableName', {sprintf(Rkonst31)});n=n+1;

        R32=horzcat(Species{5},' ','+',' ',Species{end},' ','->',' ',Species{10});RO32=horzcat('p52_processing',Species{5},num2str(n));KLO32=horzcat('kine',num2str(n));Rkonst32=horzcat('K',num2str(n));
        eval([sprintf(RO32) '=addreaction(SMV2,sprintf(R32));']); 
        eval([sprintf(KLO32) '= addkineticlaw(eval(sprintf(RO32)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst32),Parameters.NIK_mediated_p52,V,U2);']); 
        set(eval(sprintf(KLO32)),'ParameterVariableName', {sprintf(Rkonst32)});n=n+1;
%===========================================================
%Chequer Boxes of Inducibility 
%=============================
%                   
%Add Rules
%==========

RfntIkba = 'tIkba = ( Kct1 * ( ( 1 + ( W11 * ( RelAp50n / Kd11 ) ^ 3 ) + ( W12 * ( RelAp52n / Kd12 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd11 ) ^ 3 + ( RelAp52n / Kd12 ) ^ 3 ) ) ) - Kmdeg1 * tIkba';
              
RfntIkbb = 'tIkbb = ( Kct2 * ( ( 1 + ( W21 * ( RelAp50n / Kd21 ) ^ 3 ) + ( W22 * ( RelAp52n / Kd22 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd21 ) ^ 3 + ( RelAp52n / Kd22 ) ^ 3 ) ) ^ Kdel1 ) - Kmdeg3 * tIkbb';
                 
RfntIkbe = 'tIkbe = ( Kct3 * ( ( 1 + ( W31 * ( RelAp50n / Kd31 ) ^ 3 ) + ( W32 * ( RelAp52n / Kd32 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd31 ) ^ 3 + ( RelAp52n / Kd32 ) ^ 3 ) ) ^ Kdel2 ) - Kmdeg3 * tIkbe';
                
Rfntp100 = 'tp100 = ( Kct5 * ( ( 1 + ( W41 * ( RelAp50n / Kd41 ) ^ 3 ) + ( W42 * ( RelAp52n / Kd42 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd41 ) ^ 3 + ( RelAp52n / Kd42 ) ^ 3 ) ) ^ Kdel3 ) - Kmdeg5 * tp100';
                 
RfntRelA = 'tRelA = Kct7  - Kmdeg7 * tRelA';

RfntRelB = 'tRelB = ( Kct8 * ( ( 1 + ( W51 * ( RelAp50n * Kdel5 / Kd51 ) ^ 3 ) + ( W52 * ( RelAp52n * Kdel5 / Kd52 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd51 ) ^ 3 + ( RelAp52n * Kdel5 / Kd52 ) ^ 3  ) ) ) - Kmdeg8 * tRelB';

RfntR_lB = 'tR_lB = ( Kct11 * ( ( 1 + ( W61 * ( RelAp50n * Kdel5 / Kd61 ) ^ 3 ) + ( W62 * ( RelAp52n * Kdel5 / Kd62 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd61 ) ^ 3 + ( RelAp52n * Kdel5 / Kd62 ) ^ 3 ) ) ^ Kdel4 ) - Kmdeg11 * tR_lB';
                   
Rfntp50  = 'tp50  = Kct9 - Kmdeg9 * tp50';
for krow= 1:size(Parameters.Kd,1)
  for kcolumn= 1:size(Parameters.Kd,2)
      Rconst=horzcat('Kd',num2str(krow),num2str(kcolumn)); Rconst1=horzcat('W',num2str(krow),num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.Kd(krow,kcolumn),V,U);']);n=n+1;
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst1),Parameters.W(krow,kcolumn));']);n=n+1; 
  end
end
%============================================
  RlO1 = addrule(SMV2,sprintf(RfntIkba),'rate');
  RlO2 = addrule(SMV2,sprintf(RfntIkbb),'rate');
  RlO3 = addrule(SMV2,sprintf(RfntIkbe),'rate');
  RlO4 = addrule(SMV2,sprintf(Rfntp100),'rate');
  RlO5 = addrule(SMV2,sprintf(RfntRelA),'rate');
  RlO6 = addrule(SMV2,sprintf(RfntRelB),'rate');
  RlO7 = addrule(SMV2,sprintf(RfntR_lB),'rate');
  RlO8 = addrule(SMV2,sprintf(Rfntp50),'rate');
%============================================
Ca='ConstantValue';%VA='Value'

%DELAY
%=====
%kdels are 1 during euilibration
for kcolumn= 1:5
      Rconst=horzcat('Kdel',num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),1,Ca,false);']);
      n=n+1; 
end
%kdels are 0 at stimulation 0 min.
for P=1:4
 E=horzcat('Eee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(0));
                                         T=horzcat('time >= ',num2str(10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%kdels are 1 at stimulation phase once delay is over.
DELAY=[37 37 45 0]; %beta epsilon p100 RelB 
for P=1:4
 E=horzcat('Eeee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(1));
                                         T=horzcat('time >= ',num2str(DELAY(P)+10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%RelB star issue
%===============

      RconsT=horzcat('Kct',num2str(8));RconsTd=horzcat('Kmdeg',num2str(8));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsT),Parameters.monparams.Const_trnx(8),Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsTd),Parameters.monparams.mRNA_Deg(8),Ca,false);']); n=n+1; 
      EEE1=addevent(SMV2,'time >= 10000','Kct8 = 0');%EEE2=addevent(SMV2,'time >= 150000','Kmdeg8 = 0');

      RconsTT=horzcat('Kct',num2str(11));RconsTTd=horzcat('Kmdeg',num2str(11));
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTT),0,Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTTd),0,Ca,false);']); n=n+1; 
      %EEE3=addevent(SMV2,'time >= 150000','Kct11 =  3.6e-5');EEE4=addevent(SMV2,'time >= 150000','Kmdeg11 = 4.6e-3');
      E=horzcat('EEE',num2str(3));N=horzcat('Kct',num2str(11),' ','=',' ',num2str(Parameters.monparams.Const_trnx(11))); T=horzcat('time >= ',num2str(10000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);       
      E=horzcat('EEE',num2str(4));N=horzcat('Kmdeg',num2str(11),' ','=',' ',num2str(Parameters.monparams.mRNA_Deg(11))); %T=horzcat('time >= ',num2str(150000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);                                                                                                  
%Effective RelA heterodimer concentration in nucleus involved with RelB synthesis
evalues =  [1 .01];etime = [0 60];
 for Q=0:1:60
 Kdel5 = interp1(etime, evalues, Q,'pchip');E=horzcat('Eeff',num2str(Q));N=horzcat('Kdel5',' ','=',' ',num2str(Kdel5));
                                         T=horzcat('time >= ',num2str(10030+Q));%RelA effective concentration start going down after 30 mins
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end 
%SIGNAL
%======
%TNFc
%cvalues = [1 60 100 65 50 36 20 10 3];ctime = [0 5 10 15 20 25 30 360 480];
%cvalues =  [1 60 100 65 50 36 20 10 3];ctime = [0 5 10 15 20 25 30 360 480];
%======================================
%Add Parameters of Inducibility
%==============================
for Q=0:1:480
 NEMO = interp1(ctime, cvalues, Q,'pchip');E=horzcat('E',num2str(Q));N=horzcat('NEMO_ACTIVITY',' ','=',' ',num2str(NEMO));
                                         T=horzcat('time >= ',num2str(10000+Q));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end

%SIMULATE MODEL
%==============
CnfgstO = addconfigset(SMV2, 'De3');
set(CnfgstO.CompileOptions, 'UnitConversion', false);
set(CnfgstO.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO,'TimeUnits',' ');
set(CnfgstO,'SolverType','sundials');
set(CnfgstO,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO)
[t4,x4,names4]= sbiosimulate(SMV2,CnfgstO);
%==========================================================
%RelB synthesis is basal. It is induced neither by RelA, nor RelB.
%Add Model
%=========
n=1;
SMV2=sbiomodel('SMV2');
%Add species 
%===========
Species={'Ikba','Ikbb','Ikbe','Ikbd','p100','p105','RelA','RelB','p50','p52','R_lB',            ...
         'RelAp50','RelAp52','RelBp50','RelBp52','RelBp100','RelAp100','R_lBp50','R_lBp100',    ...
         'IkbaRelAp50','IkbbRelAp50','IkbeRelAp50','IkbdRelAp50',                               ...
         'IkbaRelAp52','IkbbRelAp52','IkbeRelAp52','IkbdRelAp52',                               ...
         'IkbaRelBp50','IkbbRelBp50','IkbeRelBp50','IkbdRelBp50',                               ...
         'IkbaR_lBp50','IkbbR_lBp50','IkbeR_lBp50','IkbdR_lBp50',                               ...
         'tIkba','tIkbb','tIkbe','tIkbd','tp100','tp105','tRelA','tRelB','tp50','tp52','tR_lB', ...
         'NEMO_ACTIVITY','NIK_ACTIVITY'};
Speciesn=cell(1,34);
for w=1:length(Species)
if w >= 36
continue
end
S=Species{w};S(end+1)='n'; Speciesn{w}=S;
end
I='InitialAmountUnit';U='nanomole';i='InitialAmount';C='ConstantAmount';B='BoundaryCondition';t='true';f='false';v=' ';
for w=1:length(Species)
if w<=35
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
else
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),true,sprintf(i),0,sprintf(I),sprintf(U));']);
end
     if w>=36;continue;end
     eval([sprintf(Speciesn{w}) '=addspecies(SMV2,Speciesn{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
end
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

set (NEMO_ACTIVITY, 'ConstantAmount', false);
set (NEMO_ACTIVITY, 'BoundaryCondition', true);
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'ConstantAmount', false);
set (NIK_ACTIVITY, 'BoundaryCondition', true);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

Parameters.monparams.Trnsln_rate   =  [12           12         12         0          12          0         12         12          12           0        12     ];
    
for q1=1:11
  if (q1~=8) && (q1~=11)
      %Constitutive transcription
      Rkonst=horzcat('Kct',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst),Parameters.monparams.Const_trnx(q1),V,U0);']);n=n+1; 
  
      %Contitutive mRNA degradation
      Rconst=horzcat('Kmdeg',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.monparams.mRNA_Deg(q1),V,U1);']);n=n+1;
  end
  if (q1~=4) && (q1~=10)
      %Translation Rate
      R3=horzcat(Species{q1+35},' ','->',' ',Species{q1});RO3=horzcat('Trnsln_',Species{q1},num2str(n));KLO3=horzcat('kine',num2str(n));Rkonst3=horzcat('K',num2str(n));
      eval([sprintf(RO3) '=addreaction(SMV2,sprintf(R3));']); 
      eval([sprintf(KLO3) '= addkineticlaw(eval(sprintf(RO3)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst3), Parameters.monparams.Trnsln_rate (q1),V,U1);']);
      set(eval(sprintf(KLO3)),'ParameterVariableName', {sprintf(Rkonst3)});n=n+1;
  end
      %Constitutive Degradation of Free primary Species (Cytoplasm)
      R4=horzcat(Species{q1},' ','->',' ','null');RO4=horzcat('Const_Prot_Deg_',Species{q1},num2str(n));KLO4=horzcat('kine',num2str(n));Rkonst4=horzcat('K',num2str(n));
      eval([sprintf(RO4) '=addreaction(SMV2,sprintf(R4));']); 
      eval([sprintf(KLO4) '= addkineticlaw(eval(sprintf(RO4)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst4),Parameters.monparams.CnstDegFrSpCyt(q1),V,U1);']); 
      set(eval(sprintf(KLO4)),'ParameterVariableName', {sprintf(Rkonst4)});n=n+1;
      %Constitutive Degradation of Free primary Species (Nucleus)
      R5=horzcat(Speciesn{q1},' ','->',' ','null');RO5=horzcat('Const_Prot_Deg_',Speciesn{q1},num2str(n));KLO5=horzcat('kine',num2str(n));Rkonst5=horzcat('K',num2str(n));
      eval([sprintf(RO5) '=addreaction(SMV2,sprintf(R5));']); 
      eval([sprintf(KLO5) '= addkineticlaw(eval(sprintf(RO5)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst5),Parameters.monparams.CnstDegFrSpNuc(q1),V,U1);']); 
      set(eval(sprintf(KLO5)),'ParameterVariableName', {sprintf(Rkonst5)});n=n+1;
      %NEMO Mediated Degradation of Canonical Ikbs
      if q1<=3
      R6=horzcat(Species{q1},' ','+',' ','NEMO_ACTIVITY',' ','->',' ','null');RO6=horzcat('NEMO_Mon_Deg_',Species{q1},num2str(n));KLO6=horzcat('kine',num2str(n));Rkonst6=horzcat('K',num2str(n));
      eval([sprintf(RO6) '=addreaction(SMV2,sprintf(R6));']); 
      eval([sprintf(KLO6) '= addkineticlaw(eval(sprintf(RO6)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst6),Parameters.monparams.NEMODegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO6)),'ParameterVariableName', {sprintf(Rkonst6)});n=n+1;
      end
      %NIK Mediated Degradation of Non-Canonical Ikb
      if q1==4
      R7=horzcat(Species{q1},' ','+',' ','NIK_ACTIVITY',' ','->',' ','null');RO7=horzcat('NIK_Mon_Deg_',Species{q1},num2str(n));KLO7=horzcat('kine',num2str(n));Rkonst7=horzcat('K',num2str(n));
      eval([sprintf(RO7) '=addreaction(SMV2,sprintf(R7));']); 
      eval([sprintf(KLO7) '= addkineticlaw(eval(sprintf(RO7)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst7), Parameters.monparams.NIKDegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO7)),'ParameterVariableName', {sprintf(Rkonst7)});n=n+1;
      end
      %Nuclear Import 
      if q1<=4
      R8=horzcat(Species{q1},' ','->',' ',Speciesn{q1});RO8=horzcat('Nuc_Im_Mon',Species{q1},num2str(n));KLO8=horzcat('kine',num2str(n));Rkonst8=horzcat('K',num2str(n));
      eval([sprintf(RO8) '=addreaction(SMV2,sprintf(R8));']); 
      eval([sprintf(KLO8) '= addkineticlaw(eval(sprintf(RO8)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst8),Parameters.monparams.Nuc_Imp(q1),V,U1);']); 
      set(eval(sprintf(KLO8)),'ParameterVariableName', {sprintf(Rkonst8)});n=n+1;
      end
      %Nuclear Export
      if q1<=4
      R9=horzcat(Speciesn{q1},' ','->',' ',Species{q1});RO9=horzcat('Nuc_Ex_Mon',Species{q1},num2str(n));KLO9=horzcat('kine',num2str(n));Rkonst9=horzcat('K',num2str(n));
      eval([sprintf(RO9) '=addreaction(SMV2,sprintf(R9));']); 
      eval([sprintf(KLO9) '= addkineticlaw(eval(sprintf(RO9)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst9),Parameters.monparams.Nuc_Exp(q1),V,U1);']); 
      set(eval(sprintf(KLO9)),'ParameterVariableName', {sprintf(Rkonst9)});n=n+1;
      end
end
%====================================================1.6e-3 


for q2=1:8
       %Association  of Dimers in Cytoplasm
       R10=horzcat(Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end),' ','->',' ',Species{11+q2});RO10=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO10=horzcat('kine',num2str(n));Rkonst10=horzcat('K',num2str(n));
       eval([sprintf(RO10) '=addreaction(SMV2,sprintf(R10));']); 
       eval([sprintf(KLO10) '= addkineticlaw(eval(sprintf(RO10)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst10),Parameters.dimparams.MonAssCyt(q2),V,U2);']); 
       set(eval(sprintf(KLO10)),'ParameterVariableName', {sprintf(Rkonst10)});n=n+1;
       %Association  of Dimers in Nucleus
       R11=horzcat(Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end),' ','->',' ',Speciesn{11+q2});RO11=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO11=horzcat('kine',num2str(n));Rkonst11=horzcat('K',num2str(n));
       eval([sprintf(RO11) '=addreaction(SMV2,sprintf(R11));']); 
       eval([sprintf(KLO11) '= addkineticlaw(eval(sprintf(RO11)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst11),Parameters.dimparams.MonAssNuc(q2),V,U2);']); 
       set(eval(sprintf(KLO11)),'ParameterVariableName', {sprintf(Rkonst11)});n=n+1;
       %Dissociation  of Dimers in Cytoplasm
       R12=horzcat(Species{11+q2},' ','->',' ',Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end));RO12=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO12=horzcat('kine',num2str(n));Rkonst12=horzcat('K',num2str(n));
       eval([sprintf(RO12) '=addreaction(SMV2,sprintf(R12));']); 
       eval([sprintf(KLO12) '= addkineticlaw(eval(sprintf(RO12)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst12),Parameters.dimparams.DimDissCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO12)),'ParameterVariableName', {sprintf(Rkonst12)});n=n+1;
       %Dissociation  of Dimers in Nucleus
       R13=horzcat(Speciesn{11+q2},' ','->',' ',Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end));RO13=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO13=horzcat('kine',num2str(n));Rkonst13=horzcat('K',num2str(n));
       eval([sprintf(RO13) '=addreaction(SMV2,sprintf(R13));']); 
       eval([sprintf(KLO13) '= addkineticlaw(eval(sprintf(RO13)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst13),Parameters.dimparams.DimDissNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO13)),'ParameterVariableName', {sprintf(Rkonst13)});n=n+1;
       %Nuclear Import
       R14=horzcat(Species{11+q2},' ','->',' ',Speciesn{11+q2});RO14=horzcat('Nuc_Im_Dim',Species{q2+11},num2str(n));KLO14=horzcat('kine',num2str(n));Rkonst14=horzcat('K',num2str(n));
       eval([sprintf(RO14) '=addreaction(SMV2,sprintf(R14));']); 
       eval([sprintf(KLO14) '= addkineticlaw(eval(sprintf(RO14)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst14),Parameters.dimparams.NucImp(q2),V,U1);']); 
       set(eval(sprintf(KLO14)),'ParameterVariableName', {sprintf(Rkonst14)});n=n+1;
       %Nuclear Export
       R15=horzcat(Speciesn{11+q2},' ','->',' ',Species{11+q2});RO15=horzcat('Nuc_Ex_Dim',Species{q2+11},num2str(n));KLO15=horzcat('kine',num2str(n));Rkonst15=horzcat('K',num2str(n));
       eval([sprintf(RO15) '=addreaction(SMV2,sprintf(R15));']); 
       eval([sprintf(KLO15) '= addkineticlaw(eval(sprintf(RO15)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst15),Parameters.dimparams.NucExp(q2),V,U1);']); 
       set(eval(sprintf(KLO15)),'ParameterVariableName', {sprintf(Rkonst15)});n=n+1;
       %Free Dimer Degradation in Cytoplasm
       R16=horzcat(Species{11+q2},' ','->',' ','null');RO16=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO16=horzcat('kine',num2str(n));Rkonst16=horzcat('K',num2str(n));
       eval([sprintf(RO16) '=addreaction(SMV2,sprintf(R16));']); 
       eval([sprintf(KLO16) '= addkineticlaw(eval(sprintf(RO16)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst16),Parameters.dimparams.DimDegCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO16)),'ParameterVariableName', {sprintf(Rkonst16)});n=n+1;
       %Free Dimer Degradation in Nucleus
       R17=horzcat(Speciesn{11+q2},' ','->',' ','null');RO17=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO17=horzcat('kine',num2str(n));Rkonst17=horzcat('K',num2str(n));
       eval([sprintf(RO17) '=addreaction(SMV2,sprintf(R17));']); 
       eval([sprintf(KLO17) '= addkineticlaw(eval(sprintf(RO17)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst17),Parameters.dimparams.DimDegNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO17)),'ParameterVariableName', {sprintf(Rkonst17)});n=n+1;
end
%========================================================
                                                     
for q3=1:16
        %Association of trimers in Cytoplasm
        R18=horzcat(Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end),' ','->',' ',Species{q3+19});RO18=horzcat('Trim_Ass_Cyto',Species{q3+19},num2str(n));KLO18=horzcat('kine',num2str(n));Rkonst18=horzcat('K',num2str(n));
        eval([sprintf(RO18) '=addreaction(SMV2,sprintf(R18));']); 
        eval([sprintf(KLO18) '= addkineticlaw(eval(sprintf(RO18)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst18),Parameters.trimparams.Trim_Ass_Cyt(q3),V,U2);']); 
        set(eval(sprintf(KLO18)),'ParameterVariableName', {sprintf(Rkonst18)});n=n+1;
        %Association of trimers in Nucleus
        R19=horzcat(Speciesn{q3+19}(1:4),'n',' ','+',' ',Speciesn{q3+19}(5:end),' ','->',' ',Speciesn{q3+19});RO19=horzcat('Trim_Ass_Nuc',Species{q3+19},num2str(n));KLO19=horzcat('kine',num2str(n));Rkonst19=horzcat('K',num2str(n));
        eval([sprintf(RO19) '=addreaction(SMV2,sprintf(R19));']); 
        eval([sprintf(KLO19) '= addkineticlaw(eval(sprintf(RO19)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst19),Parameters.trimparams.Trim_Ass_Nuc(q3),V,U2);']); 
        set(eval(sprintf(KLO19)),'ParameterVariableName', {sprintf(Rkonst19)});n=n+1;
        %Dissociation of trimers in Cytoplasm
        R20=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end));RO20=horzcat('Trim_Dis_Cyt',Species{q3+19},num2str(n));KLO20=horzcat('kine',num2str(n));Rkonst20=horzcat('K',num2str(n));
        eval([sprintf(RO20) '=addreaction(SMV2,sprintf(R20));']); 
        eval([sprintf(KLO20) '= addkineticlaw(eval(sprintf(RO20)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst20),Parameters.trimparams.Trim_Dis_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO20)),'ParameterVariableName', {sprintf(Rkonst20)});n=n+1;
        %Dissociation of trimers in Nucleus
        R21=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n',' ','+',' ',Species{q3+19}(5:end),'n');RO21=horzcat('Trim_Dis_Nuc',Species{q3+19},num2str(n));KLO21=horzcat('kine',num2str(n));Rkonst21=horzcat('K',num2str(n));
        eval([sprintf(RO21) '=addreaction(SMV2,sprintf(R21));']); 
        eval([sprintf(KLO21) '= addkineticlaw(eval(sprintf(RO21)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst21),Parameters.trimparams.Trim_Dis_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO21)),'ParameterVariableName', {sprintf(Rkonst21)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Cytoplasm
        R22=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(5:end));RO22=horzcat('Const_Deg_Trim_Ikb_Cyt',Species{q3+19},num2str(n));KLO22=horzcat('kine',num2str(n));Rkonst22=horzcat('K',num2str(n));
        eval([sprintf(RO22) '=addreaction(SMV2,sprintf(R22));']); 
        eval([sprintf(KLO22) '= addkineticlaw(eval(sprintf(RO22)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst22),Parameters.trimparams.Const_Trim_Deg_IkB_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO22)),'ParameterVariableName', {sprintf(Rkonst22)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Nucleus
        R23=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(5:end),'n');RO23=horzcat('Const_Deg_Trim_Ikb_Nuc',Species{q3+19},num2str(n));KLO23=horzcat('kine',num2str(n));Rkonst23=horzcat('K',num2str(n));
        eval([sprintf(RO23) '=addreaction(SMV2,sprintf(R23));']); 
        eval([sprintf(KLO23) '= addkineticlaw(eval(sprintf(RO23)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst23),Parameters.trimparams.Const_Trim_Deg_IkB_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO23)),'ParameterVariableName', {sprintf(Rkonst23)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Cytoplasm
        R24=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4));RO24=horzcat('Const_Deg_Trim_NFkb_Cyt',Species{q3+19},num2str(n));KLO24=horzcat('kine',num2str(n));Rkonst24=horzcat('K',num2str(n));
        eval([sprintf(RO24) '=addreaction(SMV2,sprintf(R24));']); 
        eval([sprintf(KLO24) '= addkineticlaw(eval(sprintf(RO24)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst24),Parameters.trimparams.Const_Trim_Deg_NFkb_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO24)),'ParameterVariableName', {sprintf(Rkonst24)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Nucleus
        R25=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n');RO25=horzcat('Const_Deg_Trim_NFkb_Nuc',Species{q3+19},num2str(n));KLO25=horzcat('kine',num2str(n));Rkonst25=horzcat('K',num2str(n));
        eval([sprintf(RO25) '=addreaction(SMV2,sprintf(R25));']); 
        eval([sprintf(KLO25) '= addkineticlaw(eval(sprintf(RO25)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst25),Parameters.trimparams.Const_Trim_Deg_NFkb_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO25)),'ParameterVariableName', {sprintf(Rkonst25)});n=n+1;
        %NEMO Mediated Degradation of Ikb in Trimer in Cytoplasm
        R26=horzcat(Species{q3+19},' ','+',' ','NEMO_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO26=horzcat('NEMO_DEg_inTrim',Species{q3+19},num2str(n));KLO26=horzcat('kine',num2str(n));Rkonst26=horzcat('K',num2str(n));
        eval([sprintf(RO26) '=addreaction(SMV2,sprintf(R26));']); 
        eval([sprintf(KLO26) '= addkineticlaw(eval(sprintf(RO26)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst26),Parameters.trimparams.NEMO_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO26)),'ParameterVariableName', {sprintf(Rkonst26)});n=n+1;
        %NIK Mediated Degradation of Ikb in Trimer in Cytoplasm
        R27=horzcat(Species{q3+19},' ','+',' ','NIK_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO27=horzcat('NIK_DEg_inTrim',Species{q3+19},num2str(n));KLO27=horzcat('kine',num2str(n));Rkonst27=horzcat('K',num2str(n));
        eval([sprintf(RO27) '=addreaction(SMV2,sprintf(R27));']); 
        eval([sprintf(KLO27) '= addkineticlaw(eval(sprintf(RO27)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst27),Parameters.trimparams.NIK_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO27)),'ParameterVariableName', {sprintf(Rkonst27)});n=n+1;
        %Nuclear Import of Trimers
        R28=horzcat(Species{q3+18},' ','->',Speciesn{q3+18});RO28=horzcat('Nuc_Imp_Trim',Species{q3+18},num2str(n));KLO28=horzcat('kine',num2str(n));Rkonst28=horzcat('K',num2str(n));
        eval([sprintf(RO28) '=addreaction(SMV2,sprintf(R28));']); 
        eval([sprintf(KLO28) '= addkineticlaw(eval(sprintf(RO28)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst28),Parameters.trimparams.Nuc_Imp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO28)),'ParameterVariableName', {sprintf(Rkonst28)});n=n+1;
        %Nuclear Export of Trimers
        R29=horzcat(Speciesn{q3+19},' ','->',Species{q3+19});RO29=horzcat('Nuc_Exp_Trip',Species{q3+19},num2str(n));KLO29=horzcat('kine',num2str(n));Rkonst29=horzcat('K',num2str(n));
        eval([sprintf(RO29) '=addreaction(SMV2,sprintf(R29));']); 
        eval([sprintf(KLO29) '= addkineticlaw(eval(sprintf(RO29)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst29),Parameters.trimparams.Nuc_Exp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO29)),'ParameterVariableName', {sprintf(Rkonst29)});n=n+1;
end
%=========================================================
%Special Reactions
   Parameters.p100_Ass            = 9e-2;     % p100 association rate to form IkBd in cytoplasm %%% 1.2E-3
   %                           K253 
   Parameters.Ikbd_Diss           = 1.2e-5;   % IkBd dissocation rate to form p100 in cytoplasm %%% 1.2E-2
   %                           K254
   Parameters.NIK_mediated_p52    = 3.5e-3;   % NIK mediated generation of P52 through processing
   %                           K255
   
        R30=horzcat('2',' ',Species{5},' ','->',Species{4});RO30=horzcat('p100_Ass',Species{5},num2str(n));KLO30=horzcat('kine',num2str(n));Rkonst30=horzcat('K',num2str(n));
        eval([sprintf(RO30) '=addreaction(SMV2,sprintf(R30));']); 
        eval([sprintf(KLO30) '= addkineticlaw(eval(sprintf(RO30)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst30),Parameters.p100_Ass,V,U2);']); 
        set(eval(sprintf(KLO30)),'ParameterVariableName', {sprintf(Rkonst30)});n=n+1;

        R31=horzcat(Species{4},' ','->',' ','2',' ',Species{5});RO31=horzcat('p100_Ass',Species{4},num2str(n));KLO31=horzcat('kine',num2str(n));Rkonst31=horzcat('K',num2str(n));
        eval([sprintf(RO31) '=addreaction(SMV2,sprintf(R31));']); 
        eval([sprintf(KLO31) '= addkineticlaw(eval(sprintf(RO31)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst31),Parameters.Ikbd_Diss,V,U1);']); 
        set(eval(sprintf(KLO31)),'ParameterVariableName', {sprintf(Rkonst31)});n=n+1;

        R32=horzcat(Species{5},' ','+',' ',Species{end},' ','->',' ',Species{10});RO32=horzcat('p52_processing',Species{5},num2str(n));KLO32=horzcat('kine',num2str(n));Rkonst32=horzcat('K',num2str(n));
        eval([sprintf(RO32) '=addreaction(SMV2,sprintf(R32));']); 
        eval([sprintf(KLO32) '= addkineticlaw(eval(sprintf(RO32)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst32),Parameters.NIK_mediated_p52,V,U2);']); 
        set(eval(sprintf(KLO32)),'ParameterVariableName', {sprintf(Rkonst32)});n=n+1;
%===========================================================
%Chequer Boxes of Inducibility 
%=============================
%                   
%Add Rules
%==========

RfntIkba = 'tIkba = ( Kct1 * ( ( 1 + ( W11 * ( RelAp50n / Kd11 ) ^ 3 ) + ( W12 * ( RelAp52n / Kd12 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd11 ) ^ 3 + ( RelAp52n / Kd12 ) ^ 3 ) ) ) - Kmdeg1 * tIkba';
              
RfntIkbb = 'tIkbb = ( Kct2 * ( ( 1 + ( W21 * ( RelAp50n / Kd21 ) ^ 3 ) + ( W22 * ( RelAp52n / Kd22 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd21 ) ^ 3 + ( RelAp52n / Kd22 ) ^ 3 ) ) ^ Kdel1 ) - Kmdeg3 * tIkbb';
                 
RfntIkbe = 'tIkbe = ( Kct3 * ( ( 1 + ( W31 * ( RelAp50n / Kd31 ) ^ 3 ) + ( W32 * ( RelAp52n / Kd32 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd31 ) ^ 3 + ( RelAp52n / Kd32 ) ^ 3 ) ) ^ Kdel2 ) - Kmdeg3 * tIkbe';
                
Rfntp100 = 'tp100 = ( Kct5 * ( ( 1 + ( W41 * ( RelAp50n / Kd41 ) ^ 3 ) + ( W42 * ( RelAp52n / Kd42 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd41 ) ^ 3 + ( RelAp52n / Kd42 ) ^ 3 ) ) ^ Kdel3 ) - Kmdeg5 * tp100';
                 
RfntRelA = 'tRelA = Kct7  - Kmdeg7 * tRelA';

RfntRelB = 'tRelB =  Kct8 - Kmdeg8 * tRelB';

RfntR_lB = 'tR_lB =  Kct11  - Kmdeg11 * tR_lB';
                   
Rfntp50  = 'tp50  = Kct9 - Kmdeg9 * tp50';
for krow= 1:size(Parameters.Kd,1)
  for kcolumn= 1:size(Parameters.Kd,2)
      Rconst=horzcat('Kd',num2str(krow),num2str(kcolumn)); Rconst1=horzcat('W',num2str(krow),num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.Kd(krow,kcolumn),V,U);']);n=n+1;
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst1),Parameters.W(krow,kcolumn));']);n=n+1; 
  end
end
%============================================
  RlO1 = addrule(SMV2,sprintf(RfntIkba),'rate');
  RlO2 = addrule(SMV2,sprintf(RfntIkbb),'rate');
  RlO3 = addrule(SMV2,sprintf(RfntIkbe),'rate');
  RlO4 = addrule(SMV2,sprintf(Rfntp100),'rate');
  RlO5 = addrule(SMV2,sprintf(RfntRelA),'rate');
  RlO6 = addrule(SMV2,sprintf(RfntRelB),'rate');
  RlO7 = addrule(SMV2,sprintf(RfntR_lB),'rate');
  RlO8 = addrule(SMV2,sprintf(Rfntp50),'rate');
%============================================
Ca='ConstantValue';%VA='Value'

%DELAY
%=====
%kdels are 1 during euilibration
for kcolumn= 1:5
      Rconst=horzcat('Kdel',num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),1,Ca,false);']);
      n=n+1; 
end
%kdels are 0 at stimulation 0 min.
for P=1:4
 E=horzcat('Eee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(0));
                                         T=horzcat('time >= ',num2str(10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%kdels are 1 at stimulation phase once delay is over.
DELAY=[37 37 45 0]; %beta epsilon p100 RelB 
for P=1:4
 E=horzcat('Eeee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(1));
                                         T=horzcat('time >= ',num2str(DELAY(P)+10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%RelB star issue
%===============

      RconsT=horzcat('Kct',num2str(8));RconsTd=horzcat('Kmdeg',num2str(8));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsT),Parameters.monparams.Const_trnx(8),Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsTd),Parameters.monparams.mRNA_Deg(8),Ca,false);']); n=n+1; 
      EEE1=addevent(SMV2,'time >= 10000','Kct8 = 0');%EEE2=addevent(SMV2,'time >= 150000','Kmdeg8 = 0');

      RconsTT=horzcat('Kct',num2str(11));RconsTTd=horzcat('Kmdeg',num2str(11));
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTT),0,Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTTd),0,Ca,false);']); n=n+1; 
      %EEE3=addevent(SMV2,'time >= 150000','Kct11 =  3.6e-5');EEE4=addevent(SMV2,'time >= 150000','Kmdeg11 = 4.6e-3');
      E=horzcat('EEE',num2str(3));N=horzcat('Kct',num2str(11),' ','=',' ',num2str(Parameters.monparams.Const_trnx(11))); T=horzcat('time >= ',num2str(10000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);       
      E=horzcat('EEE',num2str(4));N=horzcat('Kmdeg',num2str(11),' ','=',' ',num2str(Parameters.monparams.mRNA_Deg(11))); %T=horzcat('time >= ',num2str(150000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);                                                                                                  
%Effective RelA heterodimer concentration in nucleus involved with RelB synthesis
evalues =  [1 .01];etime = [0 60];
 for Q=0:1:60
 Kdel5 = interp1(etime, evalues, Q,'pchip');E=horzcat('Eeff',num2str(Q));N=horzcat('Kdel5',' ','=',' ',num2str(Kdel5));
                                         T=horzcat('time >= ',num2str(10030+Q));%RelA effective concentration start going down after 30 mins
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end 
%SIGNAL
%======
%TNFc
%cvalues = [1 60 100 65 50 36 20 10 3];ctime = [0 5 10 15 20 25 30 360 480];
%cvalues =  [1 60 100 65 50 36 20 10 3];ctime = [0 5 10 15 20 25 30 360 480];
%======================================
%Add Parameters of Inducibility
%==============================
for Q=0:1:480
 NEMO = interp1(ctime, cvalues, Q,'pchip');E=horzcat('E',num2str(Q));N=horzcat('NEMO_ACTIVITY',' ','=',' ',num2str(NEMO));
                                         T=horzcat('time >= ',num2str(10000+Q));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end

%SIMULATE MODEL
%==============
CnfgstO = addconfigset(SMV2, 'De3');
set(CnfgstO.CompileOptions, 'UnitConversion', false);
set(CnfgstO.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO,'TimeUnits',' ');
set(CnfgstO,'SolverType','sundials');
set(CnfgstO,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO)
[t5,x5,names7]= sbiosimulate(SMV2,CnfgstO);
%ADD VARIANT 1
%==============
estvarObj1 = addvariant (SMV2, 'nfkb2kn');
addcontent(estvarObj1, {'parameter','K34', 'Value',0});
CnfgstO1 = addconfigset(SMV2, 'nfkb2kn');
set(CnfgstO1.CompileOptions, 'UnitConversion', false);
set(CnfgstO1.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO1,'TimeUnits',' ');
set(CnfgstO1,'SolverType','sundials');
set(CnfgstO1,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO1,estvarObj1)
[t6,x6,names8]= sbiosimulate(SMV2,CnfgstO1,estvarObj1);




%==========================================
%==========================================
%Add Model
%=========
SMV2=sbiomodel('SMV2');
%Add species 
%===========
Species={'Ikba','Ikbb','Ikbe','Ikbd','p100','p105','RelA','RelB','p50','p52','R_lB',            ...
         'RelAp50','RelAp52','RelBp50','RelBp52','RelBp100','RelAp100','R_lBp50','R_lBp100',    ...
         'IkbaRelAp50','IkbbRelAp50','IkbeRelAp50','IkbdRelAp50',                               ...
         'IkbaRelAp52','IkbbRelAp52','IkbeRelAp52','IkbdRelAp52',                               ...
         'IkbaRelBp50','IkbbRelBp50','IkbeRelBp50','IkbdRelBp50',                               ...
         'IkbaR_lBp50','IkbbR_lBp50','IkbeR_lBp50','IkbdR_lBp50',                               ...
         'tIkba','tIkbb','tIkbe','tIkbd','tp100','tp105','tRelA','tRelB','tp50','tp52','tR_lB', ...
         'NEMO_ACTIVITY','NIK_ACTIVITY'};
Speciesn=cell(1,34);
for w=1:length(Species)
if w >= 36
continue
end
S=Species{w};S(end+1)='n'; Speciesn{w}=S;
end
I='InitialAmountUnit';U='nanomole';i='InitialAmount';C='ConstantAmount';B='BoundaryCondition';t='true';f='false';v=' ';
for w=1:length(Species)
if w<=35
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
else
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),true,sprintf(i),0,sprintf(I),sprintf(U));']);
end
     if w>=36;continue;end
     eval([sprintf(Speciesn{w}) '=addspecies(SMV2,Speciesn{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
end
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

set (NEMO_ACTIVITY, 'ConstantAmount', false);
set (NEMO_ACTIVITY, 'BoundaryCondition', true);
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'ConstantAmount', false);
set (NIK_ACTIVITY, 'BoundaryCondition', true);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

Parameters.monparams.Trnsln_rate   =  [12           12         12         0          0          0         12         12          12          0        12        ];
    
for q1=1:11
  if (q1~=8) && (q1~=11)
      %Constitutive transcription
      Rkonst=horzcat('Kct',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst),Parameters.monparams.Const_trnx(q1),V,U0);']);n=n+1; 
  
      %Contitutive mRNA degradation
      Rconst=horzcat('Kmdeg',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.monparams.mRNA_Deg(q1),V,U1);']);n=n+1;
  end
  if (q1~=4) && (q1~=10)
      %Translation Rate
      R3=horzcat(Species{q1+35},' ','->',' ',Species{q1});RO3=horzcat('Trnsln_',Species{q1},num2str(n));KLO3=horzcat('kine',num2str(n));Rkonst3=horzcat('K',num2str(n));
      eval([sprintf(RO3) '=addreaction(SMV2,sprintf(R3));']); 
      eval([sprintf(KLO3) '= addkineticlaw(eval(sprintf(RO3)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst3), Parameters.monparams.Trnsln_rate (q1),V,U1);']);
      set(eval(sprintf(KLO3)),'ParameterVariableName', {sprintf(Rkonst3)});n=n+1;
  end
      %Constitutive Degradation of Free primary Species (Cytoplasm)
      R4=horzcat(Species{q1},' ','->',' ','null');RO4=horzcat('Const_Prot_Deg_',Species{q1},num2str(n));KLO4=horzcat('kine',num2str(n));Rkonst4=horzcat('K',num2str(n));
      eval([sprintf(RO4) '=addreaction(SMV2,sprintf(R4));']); 
      eval([sprintf(KLO4) '= addkineticlaw(eval(sprintf(RO4)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst4),Parameters.monparams.CnstDegFrSpCyt(q1),V,U1);']); 
      set(eval(sprintf(KLO4)),'ParameterVariableName', {sprintf(Rkonst4)});n=n+1;
      %Constitutive Degradation of Free primary Species (Nucleus)
      R5=horzcat(Speciesn{q1},' ','->',' ','null');RO5=horzcat('Const_Prot_Deg_',Speciesn{q1},num2str(n));KLO5=horzcat('kine',num2str(n));Rkonst5=horzcat('K',num2str(n));
      eval([sprintf(RO5) '=addreaction(SMV2,sprintf(R5));']); 
      eval([sprintf(KLO5) '= addkineticlaw(eval(sprintf(RO5)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst5),Parameters.monparams.CnstDegFrSpNuc(q1),V,U1);']); 
      set(eval(sprintf(KLO5)),'ParameterVariableName', {sprintf(Rkonst5)});n=n+1;
      %NEMO Mediated Degradation of Canonical Ikbs
      if q1<=3
      R6=horzcat(Species{q1},' ','+',' ','NEMO_ACTIVITY',' ','->',' ','null');RO6=horzcat('NEMO_Mon_Deg_',Species{q1},num2str(n));KLO6=horzcat('kine',num2str(n));Rkonst6=horzcat('K',num2str(n));
      eval([sprintf(RO6) '=addreaction(SMV2,sprintf(R6));']); 
      eval([sprintf(KLO6) '= addkineticlaw(eval(sprintf(RO6)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst6),Parameters.monparams.NEMODegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO6)),'ParameterVariableName', {sprintf(Rkonst6)});n=n+1;
      end
      %NIK Mediated Degradation of Non-Canonical Ikb
      if q1==4
      R7=horzcat(Species{q1},' ','+',' ','NIK_ACTIVITY',' ','->',' ','null');RO7=horzcat('NIK_Mon_Deg_',Species{q1},num2str(n));KLO7=horzcat('kine',num2str(n));Rkonst7=horzcat('K',num2str(n));
      eval([sprintf(RO7) '=addreaction(SMV2,sprintf(R7));']); 
      eval([sprintf(KLO7) '= addkineticlaw(eval(sprintf(RO7)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst7), Parameters.monparams.NIKDegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO7)),'ParameterVariableName', {sprintf(Rkonst7)});n=n+1;
      end
      %Nuclear Import 
      if q1<=4
      R8=horzcat(Species{q1},' ','->',' ',Speciesn{q1});RO8=horzcat('Nuc_Im_Mon',Species{q1},num2str(n));KLO8=horzcat('kine',num2str(n));Rkonst8=horzcat('K',num2str(n));
      eval([sprintf(RO8) '=addreaction(SMV2,sprintf(R8));']); 
      eval([sprintf(KLO8) '= addkineticlaw(eval(sprintf(RO8)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst8),Parameters.monparams.Nuc_Imp(q1),V,U1);']); 
      set(eval(sprintf(KLO8)),'ParameterVariableName', {sprintf(Rkonst8)});n=n+1;
      end
      %Nuclear Export
      if q1<=4
      R9=horzcat(Speciesn{q1},' ','->',' ',Species{q1});RO9=horzcat('Nuc_Ex_Mon',Species{q1},num2str(n));KLO9=horzcat('kine',num2str(n));Rkonst9=horzcat('K',num2str(n));
      eval([sprintf(RO9) '=addreaction(SMV2,sprintf(R9));']); 
      eval([sprintf(KLO9) '= addkineticlaw(eval(sprintf(RO9)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst9),Parameters.monparams.Nuc_Exp(q1),V,U1);']); 
      set(eval(sprintf(KLO9)),'ParameterVariableName', {sprintf(Rkonst9)});n=n+1;
      end
end
%====================================================1.6e-3 


for q2=1:8
       %Association  of Dimers in Cytoplasm
       R10=horzcat(Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end),' ','->',' ',Species{11+q2});RO10=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO10=horzcat('kine',num2str(n));Rkonst10=horzcat('K',num2str(n));
       eval([sprintf(RO10) '=addreaction(SMV2,sprintf(R10));']); 
       eval([sprintf(KLO10) '= addkineticlaw(eval(sprintf(RO10)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst10),Parameters.dimparams.MonAssCyt(q2),V,U2);']); 
       set(eval(sprintf(KLO10)),'ParameterVariableName', {sprintf(Rkonst10)});n=n+1;
       %Association  of Dimers in Nucleus
       R11=horzcat(Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end),' ','->',' ',Speciesn{11+q2});RO11=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO11=horzcat('kine',num2str(n));Rkonst11=horzcat('K',num2str(n));
       eval([sprintf(RO11) '=addreaction(SMV2,sprintf(R11));']); 
       eval([sprintf(KLO11) '= addkineticlaw(eval(sprintf(RO11)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst11),Parameters.dimparams.MonAssNuc(q2),V,U2);']); 
       set(eval(sprintf(KLO11)),'ParameterVariableName', {sprintf(Rkonst11)});n=n+1;
       %Dissociation  of Dimers in Cytoplasm
       R12=horzcat(Species{11+q2},' ','->',' ',Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end));RO12=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO12=horzcat('kine',num2str(n));Rkonst12=horzcat('K',num2str(n));
       eval([sprintf(RO12) '=addreaction(SMV2,sprintf(R12));']); 
       eval([sprintf(KLO12) '= addkineticlaw(eval(sprintf(RO12)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst12),Parameters.dimparams.DimDissCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO12)),'ParameterVariableName', {sprintf(Rkonst12)});n=n+1;
       %Dissociation  of Dimers in Nucleus
       R13=horzcat(Speciesn{11+q2},' ','->',' ',Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end));RO13=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO13=horzcat('kine',num2str(n));Rkonst13=horzcat('K',num2str(n));
       eval([sprintf(RO13) '=addreaction(SMV2,sprintf(R13));']); 
       eval([sprintf(KLO13) '= addkineticlaw(eval(sprintf(RO13)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst13),Parameters.dimparams.DimDissNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO13)),'ParameterVariableName', {sprintf(Rkonst13)});n=n+1;
       %Nuclear Import
       R14=horzcat(Species{11+q2},' ','->',' ',Speciesn{11+q2});RO14=horzcat('Nuc_Im_Dim',Species{q2+11},num2str(n));KLO14=horzcat('kine',num2str(n));Rkonst14=horzcat('K',num2str(n));
       eval([sprintf(RO14) '=addreaction(SMV2,sprintf(R14));']); 
       eval([sprintf(KLO14) '= addkineticlaw(eval(sprintf(RO14)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst14),Parameters.dimparams.NucImp(q2),V,U1);']); 
       set(eval(sprintf(KLO14)),'ParameterVariableName', {sprintf(Rkonst14)});n=n+1;
       %Nuclear Export
       R15=horzcat(Speciesn{11+q2},' ','->',' ',Species{11+q2});RO15=horzcat('Nuc_Ex_Dim',Species{q2+11},num2str(n));KLO15=horzcat('kine',num2str(n));Rkonst15=horzcat('K',num2str(n));
       eval([sprintf(RO15) '=addreaction(SMV2,sprintf(R15));']); 
       eval([sprintf(KLO15) '= addkineticlaw(eval(sprintf(RO15)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst15),Parameters.dimparams.NucExp(q2),V,U1);']); 
       set(eval(sprintf(KLO15)),'ParameterVariableName', {sprintf(Rkonst15)});n=n+1;
       %Free Dimer Degradation in Cytoplasm
       R16=horzcat(Species{11+q2},' ','->',' ','null');RO16=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO16=horzcat('kine',num2str(n));Rkonst16=horzcat('K',num2str(n));
       eval([sprintf(RO16) '=addreaction(SMV2,sprintf(R16));']); 
       eval([sprintf(KLO16) '= addkineticlaw(eval(sprintf(RO16)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst16),Parameters.dimparams.DimDegCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO16)),'ParameterVariableName', {sprintf(Rkonst16)});n=n+1;
       %Free Dimer Degradation in Nucleus
       R17=horzcat(Speciesn{11+q2},' ','->',' ','null');RO17=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO17=horzcat('kine',num2str(n));Rkonst17=horzcat('K',num2str(n));
       eval([sprintf(RO17) '=addreaction(SMV2,sprintf(R17));']); 
       eval([sprintf(KLO17) '= addkineticlaw(eval(sprintf(RO17)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst17),Parameters.dimparams.DimDegNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO17)),'ParameterVariableName', {sprintf(Rkonst17)});n=n+1;
end
%========================================================
                                                     
for q3=1:16
        %Association of trimers in Cytoplasm
        R18=horzcat(Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end),' ','->',' ',Species{q3+19});RO18=horzcat('Trim_Ass_Cyto',Species{q3+19},num2str(n));KLO18=horzcat('kine',num2str(n));Rkonst18=horzcat('K',num2str(n));
        eval([sprintf(RO18) '=addreaction(SMV2,sprintf(R18));']); 
        eval([sprintf(KLO18) '= addkineticlaw(eval(sprintf(RO18)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst18),Parameters.trimparams.Trim_Ass_Cyt(q3),V,U2);']); 
        set(eval(sprintf(KLO18)),'ParameterVariableName', {sprintf(Rkonst18)});n=n+1;
        %Association of trimers in Nucleus
        R19=horzcat(Speciesn{q3+19}(1:4),'n',' ','+',' ',Speciesn{q3+19}(5:end),' ','->',' ',Speciesn{q3+19});RO19=horzcat('Trim_Ass_Nuc',Species{q3+19},num2str(n));KLO19=horzcat('kine',num2str(n));Rkonst19=horzcat('K',num2str(n));
        eval([sprintf(RO19) '=addreaction(SMV2,sprintf(R19));']); 
        eval([sprintf(KLO19) '= addkineticlaw(eval(sprintf(RO19)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst19),Parameters.trimparams.Trim_Ass_Nuc(q3),V,U2);']); 
        set(eval(sprintf(KLO19)),'ParameterVariableName', {sprintf(Rkonst19)});n=n+1;
        %Dissociation of trimers in Cytoplasm
        R20=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end));RO20=horzcat('Trim_Dis_Cyt',Species{q3+19},num2str(n));KLO20=horzcat('kine',num2str(n));Rkonst20=horzcat('K',num2str(n));
        eval([sprintf(RO20) '=addreaction(SMV2,sprintf(R20));']); 
        eval([sprintf(KLO20) '= addkineticlaw(eval(sprintf(RO20)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst20),Parameters.trimparams.Trim_Dis_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO20)),'ParameterVariableName', {sprintf(Rkonst20)});n=n+1;
        %Dissociation of trimers in Nucleus
        R21=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n',' ','+',' ',Species{q3+19}(5:end),'n');RO21=horzcat('Trim_Dis_Nuc',Species{q3+19},num2str(n));KLO21=horzcat('kine',num2str(n));Rkonst21=horzcat('K',num2str(n));
        eval([sprintf(RO21) '=addreaction(SMV2,sprintf(R21));']); 
        eval([sprintf(KLO21) '= addkineticlaw(eval(sprintf(RO21)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst21),Parameters.trimparams.Trim_Dis_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO21)),'ParameterVariableName', {sprintf(Rkonst21)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Cytoplasm
        R22=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(5:end));RO22=horzcat('Const_Deg_Trim_Ikb_Cyt',Species{q3+19},num2str(n));KLO22=horzcat('kine',num2str(n));Rkonst22=horzcat('K',num2str(n));
        eval([sprintf(RO22) '=addreaction(SMV2,sprintf(R22));']); 
        eval([sprintf(KLO22) '= addkineticlaw(eval(sprintf(RO22)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst22),Parameters.trimparams.Const_Trim_Deg_IkB_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO22)),'ParameterVariableName', {sprintf(Rkonst22)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Nucleus
        R23=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(5:end),'n');RO23=horzcat('Const_Deg_Trim_Ikb_Nuc',Species{q3+19},num2str(n));KLO23=horzcat('kine',num2str(n));Rkonst23=horzcat('K',num2str(n));
        eval([sprintf(RO23) '=addreaction(SMV2,sprintf(R23));']); 
        eval([sprintf(KLO23) '= addkineticlaw(eval(sprintf(RO23)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst23),Parameters.trimparams.Const_Trim_Deg_IkB_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO23)),'ParameterVariableName', {sprintf(Rkonst23)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Cytoplasm
        R24=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4));RO24=horzcat('Const_Deg_Trim_NFkb_Cyt',Species{q3+19},num2str(n));KLO24=horzcat('kine',num2str(n));Rkonst24=horzcat('K',num2str(n));
        eval([sprintf(RO24) '=addreaction(SMV2,sprintf(R24));']); 
        eval([sprintf(KLO24) '= addkineticlaw(eval(sprintf(RO24)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst24),Parameters.trimparams.Const_Trim_Deg_NFkb_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO24)),'ParameterVariableName', {sprintf(Rkonst24)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Nucleus
        R25=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n');RO25=horzcat('Const_Deg_Trim_NFkb_Nuc',Species{q3+19},num2str(n));KLO25=horzcat('kine',num2str(n));Rkonst25=horzcat('K',num2str(n));
        eval([sprintf(RO25) '=addreaction(SMV2,sprintf(R25));']); 
        eval([sprintf(KLO25) '= addkineticlaw(eval(sprintf(RO25)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst25),Parameters.trimparams.Const_Trim_Deg_NFkb_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO25)),'ParameterVariableName', {sprintf(Rkonst25)});n=n+1;
        %NEMO Mediated Degradation of Ikb in Trimer in Cytoplasm
        R26=horzcat(Species{q3+19},' ','+',' ','NEMO_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO26=horzcat('NEMO_DEg_inTrim',Species{q3+19},num2str(n));KLO26=horzcat('kine',num2str(n));Rkonst26=horzcat('K',num2str(n));
        eval([sprintf(RO26) '=addreaction(SMV2,sprintf(R26));']); 
        eval([sprintf(KLO26) '= addkineticlaw(eval(sprintf(RO26)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst26),Parameters.trimparams.NEMO_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO26)),'ParameterVariableName', {sprintf(Rkonst26)});n=n+1;
        %NIK Mediated Degradation of Ikb in Trimer in Cytoplasm
        R27=horzcat(Species{q3+19},' ','+',' ','NIK_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO27=horzcat('NIK_DEg_inTrim',Species{q3+19},num2str(n));KLO27=horzcat('kine',num2str(n));Rkonst27=horzcat('K',num2str(n));
        eval([sprintf(RO27) '=addreaction(SMV2,sprintf(R27));']); 
        eval([sprintf(KLO27) '= addkineticlaw(eval(sprintf(RO27)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst27),Parameters.trimparams.NIK_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO27)),'ParameterVariableName', {sprintf(Rkonst27)});n=n+1;
        %Nuclear Import of Trimers
        R28=horzcat(Species{q3+18},' ','->',Speciesn{q3+18});RO28=horzcat('Nuc_Imp_Trim',Species{q3+18},num2str(n));KLO28=horzcat('kine',num2str(n));Rkonst28=horzcat('K',num2str(n));
        eval([sprintf(RO28) '=addreaction(SMV2,sprintf(R28));']); 
        eval([sprintf(KLO28) '= addkineticlaw(eval(sprintf(RO28)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst28),Parameters.trimparams.Nuc_Imp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO28)),'ParameterVariableName', {sprintf(Rkonst28)});n=n+1;
        %Nuclear Export of Trimers
        R29=horzcat(Speciesn{q3+19},' ','->',Species{q3+19});RO29=horzcat('Nuc_Exp_Trip',Species{q3+19},num2str(n));KLO29=horzcat('kine',num2str(n));Rkonst29=horzcat('K',num2str(n));
        eval([sprintf(RO29) '=addreaction(SMV2,sprintf(R29));']); 
        eval([sprintf(KLO29) '= addkineticlaw(eval(sprintf(RO29)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst29),Parameters.trimparams.Nuc_Exp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO29)),'ParameterVariableName', {sprintf(Rkonst29)});n=n+1;
end
%=========================================================
%Special Reactions
   Parameters.p100_Ass            = 9e-2;     % p100 association rate to form IkBd in cytoplasm %%% 1.2E-3
   %                           K253 
   Parameters.Ikbd_Diss           = 1.2e-5;   % IkBd dissocation rate to form p100 in cytoplasm %%% 1.2E-2
   %                           K254
   Parameters.NIK_mediated_p52    = 3.5e-3;   % NIK mediated generation of P52 through processing
   %                           K255
   
        R30=horzcat('2',' ',Species{5},' ','->',Species{4});RO30=horzcat('p100_Ass',Species{5},num2str(n));KLO30=horzcat('kine',num2str(n));Rkonst30=horzcat('K',num2str(n));
        eval([sprintf(RO30) '=addreaction(SMV2,sprintf(R30));']); 
        eval([sprintf(KLO30) '= addkineticlaw(eval(sprintf(RO30)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst30),Parameters.p100_Ass,V,U2);']); 
        set(eval(sprintf(KLO30)),'ParameterVariableName', {sprintf(Rkonst30)});n=n+1;

        R31=horzcat(Species{4},' ','->',' ','2',' ',Species{5});RO31=horzcat('p100_Ass',Species{4},num2str(n));KLO31=horzcat('kine',num2str(n));Rkonst31=horzcat('K',num2str(n));
        eval([sprintf(RO31) '=addreaction(SMV2,sprintf(R31));']); 
        eval([sprintf(KLO31) '= addkineticlaw(eval(sprintf(RO31)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst31),Parameters.Ikbd_Diss,V,U1);']); 
        set(eval(sprintf(KLO31)),'ParameterVariableName', {sprintf(Rkonst31)});n=n+1;

        R32=horzcat(Species{5},' ','+',' ',Species{end},' ','->',' ',Species{10});RO32=horzcat('p52_processing',Species{5},num2str(n));KLO32=horzcat('kine',num2str(n));Rkonst32=horzcat('K',num2str(n));
        eval([sprintf(RO32) '=addreaction(SMV2,sprintf(R32));']); 
        eval([sprintf(KLO32) '= addkineticlaw(eval(sprintf(RO32)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst32),Parameters.NIK_mediated_p52,V,U2);']); 
        set(eval(sprintf(KLO32)),'ParameterVariableName', {sprintf(Rkonst32)});n=n+1;
%===========================================================
%Chequer Boxes of Inducibility 
%=============================
%                   
%Add Rules
%==========

RfntIkba = 'tIkba = ( Kct1 * ( ( 1 + ( W11 * ( RelAp50n / Kd11 ) ^ 3 ) + ( W12 * ( RelAp52n / Kd12 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd11 ) ^ 3 + ( RelAp52n / Kd12 ) ^ 3 ) ) ) - Kmdeg1 * tIkba';
              
RfntIkbb = 'tIkbb = ( Kct2 * ( ( 1 + ( W21 * ( RelAp50n / Kd21 ) ^ 3 ) + ( W22 * ( RelAp52n / Kd22 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd21 ) ^ 3 + ( RelAp52n / Kd22 ) ^ 3 ) ) ^ Kdel1 ) - Kmdeg3 * tIkbb';
                 
RfntIkbe = 'tIkbe = ( Kct3 * ( ( 1 + ( W31 * ( RelAp50n / Kd31 ) ^ 3 ) + ( W32 * ( RelAp52n / Kd32 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd31 ) ^ 3 + ( RelAp52n / Kd32 ) ^ 3 ) ) ^ Kdel2 ) - Kmdeg3 * tIkbe';
                
Rfntp100 = 'tp100 = ( Kct5 * ( ( 1 + ( W41 * ( RelAp50n / Kd41 ) ^ 3 ) + ( W42 * ( RelAp52n / Kd42 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd41 ) ^ 3 + ( RelAp52n / Kd42 ) ^ 3 ) ) ^ Kdel3 ) - Kmdeg5 * tp100';
                 
RfntRelA = 'tRelA = Kct7  - Kmdeg7 * tRelA';

RfntRelB = 'tRelB = ( Kct8 * ( ( 1 + ( W51 * ( RelAp50n * Kdel5 / Kd51 ) ^ 3 ) + ( W52 * ( RelAp52n * Kdel5 / Kd52 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd51 ) ^ 3 + ( RelAp52n * Kdel5 / Kd52 ) ^ 3  ) ) ) - Kmdeg8 * tRelB';

RfntR_lB = 'tR_lB = ( Kct11 * ( ( 1 + ( W61 * ( RelAp50n * Kdel5 / Kd61 ) ^ 3 ) + ( W62 * ( RelAp52n * Kdel5 / Kd62 ) ^ 3 ) ) / ( 1 + ( RelAp50n * Kdel5 / Kd61 ) ^ 3 + ( RelAp52n * Kdel5 / Kd62 ) ^ 3 ) ) ^ Kdel4 ) - Kmdeg11 * tR_lB';
                   
Rfntp50  = 'tp50  = Kct9 - Kmdeg9 * tp50';
for krow= 1:size(Parameters.Kd,1)
  for kcolumn= 1:size(Parameters.Kd,2)
      Rconst=horzcat('Kd',num2str(krow),num2str(kcolumn)); Rconst1=horzcat('W',num2str(krow),num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.Kd(krow,kcolumn),V,U);']);n=n+1;
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst1),Parameters.W(krow,kcolumn));']);n=n+1; 
  end
end
%============================================
  RlO1 = addrule(SMV2,sprintf(RfntIkba),'rate');
  RlO2 = addrule(SMV2,sprintf(RfntIkbb),'rate');
  RlO3 = addrule(SMV2,sprintf(RfntIkbe),'rate');
  RlO4 = addrule(SMV2,sprintf(Rfntp100),'rate');
  RlO5 = addrule(SMV2,sprintf(RfntRelA),'rate');
  RlO6 = addrule(SMV2,sprintf(RfntRelB),'rate');
  RlO7 = addrule(SMV2,sprintf(RfntR_lB),'rate');
  RlO8 = addrule(SMV2,sprintf(Rfntp50),'rate');
%============================================
Ca='ConstantValue';%VA='Value'

%DELAY
%=====
%kdels are 1 during euilibration
for kcolumn= 1:5
      Rconst=horzcat('Kdel',num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),1,Ca,false);']);
      n=n+1; 
end
%kdels are 0 at stimulation 0 min.
for P=1:4
 E=horzcat('Eee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(0));
                                         T=horzcat('time >= ',num2str(10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%kdels are 1 at stimulation phase once delay is over.
DELAY=[37 37 45 0]; %beta epsilon p100 RelB 
for P=1:4
 E=horzcat('Eeee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(1));
                                         T=horzcat('time >= ',num2str(DELAY(P)+10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%RelB star issue
%===============

      RconsT=horzcat('Kct',num2str(8));RconsTd=horzcat('Kmdeg',num2str(8));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsT),Parameters.monparams.Const_trnx(8),Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsTd),Parameters.monparams.mRNA_Deg(8),Ca,false);']); n=n+1; 
      EEE1=addevent(SMV2,'time >= 10000','Kct8 = 0');%EEE2=addevent(SMV2,'time >= 150000','Kmdeg8 = 0');

      RconsTT=horzcat('Kct',num2str(11));RconsTTd=horzcat('Kmdeg',num2str(11));
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTT),0,Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTTd),0,Ca,false);']); n=n+1; 
      %EEE3=addevent(SMV2,'time >= 150000','Kct11 =  3.6e-5');EEE4=addevent(SMV2,'time >= 150000','Kmdeg11 = 4.6e-3');
      E=horzcat('EEE',num2str(3));N=horzcat('Kct',num2str(11),' ','=',' ',num2str(Parameters.monparams.Const_trnx(11))); T=horzcat('time >= ',num2str(10000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);       
      E=horzcat('EEE',num2str(4));N=horzcat('Kmdeg',num2str(11),' ','=',' ',num2str(Parameters.monparams.mRNA_Deg(11))); %T=horzcat('time >= ',num2str(150000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);                                                                                                  
%Effective RelA heterodimer concentration in nucleus involved with RelB synthesis
evalues =  [1 .01];etime = [0 60];
 for Q=0:1:60
 Kdel5 = interp1(etime, evalues, Q,'pchip');E=horzcat('Eeff',num2str(Q));N=horzcat('Kdel5',' ','=',' ',num2str(Kdel5));
                                         T=horzcat('time >= ',num2str(10030+Q));%RelA effective concentration start going down after 30 mins
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end 

for Q=0:1:480
 NEMO = interp1(ptime, pvalues, Q,'pchip');E=horzcat('E',num2str(Q));N=horzcat('NEMO_ACTIVITY',' ','=',' ',num2str(NEMO));
                                         T=horzcat('time >= ',num2str(10000+Q));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end

%SIMULATE MODEL
%==============
CnfgstO = addconfigset(SMV2, 'De3');
set(CnfgstO.CompileOptions, 'UnitConversion', false);
set(CnfgstO.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO,'TimeUnits',' ');
set(CnfgstO,'SolverType','sundials');
set(CnfgstO,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO)
[t10,x10,names6]= sbiosimulate(SMV2,CnfgstO);
%==========================================================
%RelB synthesis is basal. It is induced neither by RelA, nor RelB.
%Add Model
%=========
n=1;
SMV2=sbiomodel('SMV2');
%Add species 
%===========
Species={'Ikba','Ikbb','Ikbe','Ikbd','p100','p105','RelA','RelB','p50','p52','R_lB',            ...
         'RelAp50','RelAp52','RelBp50','RelBp52','RelBp100','RelAp100','R_lBp50','R_lBp100',    ...
         'IkbaRelAp50','IkbbRelAp50','IkbeRelAp50','IkbdRelAp50',                               ...
         'IkbaRelAp52','IkbbRelAp52','IkbeRelAp52','IkbdRelAp52',                               ...
         'IkbaRelBp50','IkbbRelBp50','IkbeRelBp50','IkbdRelBp50',                               ...
         'IkbaR_lBp50','IkbbR_lBp50','IkbeR_lBp50','IkbdR_lBp50',                               ...
         'tIkba','tIkbb','tIkbe','tIkbd','tp100','tp105','tRelA','tRelB','tp50','tp52','tR_lB', ...
         'NEMO_ACTIVITY','NIK_ACTIVITY'};
Speciesn=cell(1,34);
for w=1:length(Species)
if w >= 36
continue
end
S=Species{w};S(end+1)='n'; Speciesn{w}=S;
end
I='InitialAmountUnit';U='nanomole';i='InitialAmount';C='ConstantAmount';B='BoundaryCondition';t='true';f='false';v=' ';
for w=1:length(Species)
if w<=35
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
else
     eval([sprintf(Species{w}) '=addspecies(SMV2,Species{w},sprintf(C),false,sprintf(B),true,sprintf(i),0,sprintf(I),sprintf(U));']);
end
     if w>=36;continue;end
     eval([sprintf(Speciesn{w}) '=addspecies(SMV2,Speciesn{w},sprintf(C),false,sprintf(B),false,sprintf(i),0,sprintf(I),sprintf(U));']);
end
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

set (NEMO_ACTIVITY, 'ConstantAmount', false);
set (NEMO_ACTIVITY, 'BoundaryCondition', true);
set (NEMO_ACTIVITY, 'InitialAmount',1.00);
set (NIK_ACTIVITY, 'ConstantAmount', false);
set (NIK_ACTIVITY, 'BoundaryCondition', true);
set (NIK_ACTIVITY, 'InitialAmount',1.00);

Parameters.monparams.Trnsln_rate   =  [12           12         12         0          12          0         12         12          12           0        12     ];
    
for q1=1:11
  if (q1~=8) && (q1~=11)
      %Constitutive transcription
      Rkonst=horzcat('Kct',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst),Parameters.monparams.Const_trnx(q1),V,U0);']);n=n+1; 
  
      %Contitutive mRNA degradation
      Rconst=horzcat('Kmdeg',num2str(q1));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.monparams.mRNA_Deg(q1),V,U1);']);n=n+1;
  end
  if (q1~=4) && (q1~=10)
      %Translation Rate
      R3=horzcat(Species{q1+35},' ','->',' ',Species{q1});RO3=horzcat('Trnsln_',Species{q1},num2str(n));KLO3=horzcat('kine',num2str(n));Rkonst3=horzcat('K',num2str(n));
      eval([sprintf(RO3) '=addreaction(SMV2,sprintf(R3));']); 
      eval([sprintf(KLO3) '= addkineticlaw(eval(sprintf(RO3)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst3), Parameters.monparams.Trnsln_rate (q1),V,U1);']);
      set(eval(sprintf(KLO3)),'ParameterVariableName', {sprintf(Rkonst3)});n=n+1;
  end
      %Constitutive Degradation of Free primary Species (Cytoplasm)
      R4=horzcat(Species{q1},' ','->',' ','null');RO4=horzcat('Const_Prot_Deg_',Species{q1},num2str(n));KLO4=horzcat('kine',num2str(n));Rkonst4=horzcat('K',num2str(n));
      eval([sprintf(RO4) '=addreaction(SMV2,sprintf(R4));']); 
      eval([sprintf(KLO4) '= addkineticlaw(eval(sprintf(RO4)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst4),Parameters.monparams.CnstDegFrSpCyt(q1),V,U1);']); 
      set(eval(sprintf(KLO4)),'ParameterVariableName', {sprintf(Rkonst4)});n=n+1;
      %Constitutive Degradation of Free primary Species (Nucleus)
      R5=horzcat(Speciesn{q1},' ','->',' ','null');RO5=horzcat('Const_Prot_Deg_',Speciesn{q1},num2str(n));KLO5=horzcat('kine',num2str(n));Rkonst5=horzcat('K',num2str(n));
      eval([sprintf(RO5) '=addreaction(SMV2,sprintf(R5));']); 
      eval([sprintf(KLO5) '= addkineticlaw(eval(sprintf(RO5)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst5),Parameters.monparams.CnstDegFrSpNuc(q1),V,U1);']); 
      set(eval(sprintf(KLO5)),'ParameterVariableName', {sprintf(Rkonst5)});n=n+1;
      %NEMO Mediated Degradation of Canonical Ikbs
      if q1<=3
      R6=horzcat(Species{q1},' ','+',' ','NEMO_ACTIVITY',' ','->',' ','null');RO6=horzcat('NEMO_Mon_Deg_',Species{q1},num2str(n));KLO6=horzcat('kine',num2str(n));Rkonst6=horzcat('K',num2str(n));
      eval([sprintf(RO6) '=addreaction(SMV2,sprintf(R6));']); 
      eval([sprintf(KLO6) '= addkineticlaw(eval(sprintf(RO6)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst6),Parameters.monparams.NEMODegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO6)),'ParameterVariableName', {sprintf(Rkonst6)});n=n+1;
      end
      %NIK Mediated Degradation of Non-Canonical Ikb
      if q1==4
      R7=horzcat(Species{q1},' ','+',' ','NIK_ACTIVITY',' ','->',' ','null');RO7=horzcat('NIK_Mon_Deg_',Species{q1},num2str(n));KLO7=horzcat('kine',num2str(n));Rkonst7=horzcat('K',num2str(n));
      eval([sprintf(RO7) '=addreaction(SMV2,sprintf(R7));']); 
      eval([sprintf(KLO7) '= addkineticlaw(eval(sprintf(RO7)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst7), Parameters.monparams.NIKDegFrSp(q1),V,U2);']); 
      set(eval(sprintf(KLO7)),'ParameterVariableName', {sprintf(Rkonst7)});n=n+1;
      end
      %Nuclear Import 
      if q1<=4
      R8=horzcat(Species{q1},' ','->',' ',Speciesn{q1});RO8=horzcat('Nuc_Im_Mon',Species{q1},num2str(n));KLO8=horzcat('kine',num2str(n));Rkonst8=horzcat('K',num2str(n));
      eval([sprintf(RO8) '=addreaction(SMV2,sprintf(R8));']); 
      eval([sprintf(KLO8) '= addkineticlaw(eval(sprintf(RO8)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst8),Parameters.monparams.Nuc_Imp(q1),V,U1);']); 
      set(eval(sprintf(KLO8)),'ParameterVariableName', {sprintf(Rkonst8)});n=n+1;
      end
      %Nuclear Export
      if q1<=4
      R9=horzcat(Speciesn{q1},' ','->',' ',Species{q1});RO9=horzcat('Nuc_Ex_Mon',Species{q1},num2str(n));KLO9=horzcat('kine',num2str(n));Rkonst9=horzcat('K',num2str(n));
      eval([sprintf(RO9) '=addreaction(SMV2,sprintf(R9));']); 
      eval([sprintf(KLO9) '= addkineticlaw(eval(sprintf(RO9)),sprintf(M));']);
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst9),Parameters.monparams.Nuc_Exp(q1),V,U1);']); 
      set(eval(sprintf(KLO9)),'ParameterVariableName', {sprintf(Rkonst9)});n=n+1;
      end
end
%====================================================1.6e-3 


for q2=1:8
       %Association  of Dimers in Cytoplasm
       R10=horzcat(Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end),' ','->',' ',Species{11+q2});RO10=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO10=horzcat('kine',num2str(n));Rkonst10=horzcat('K',num2str(n));
       eval([sprintf(RO10) '=addreaction(SMV2,sprintf(R10));']); 
       eval([sprintf(KLO10) '= addkineticlaw(eval(sprintf(RO10)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst10),Parameters.dimparams.MonAssCyt(q2),V,U2);']); 
       set(eval(sprintf(KLO10)),'ParameterVariableName', {sprintf(Rkonst10)});n=n+1;
       %Association  of Dimers in Nucleus
       R11=horzcat(Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end),' ','->',' ',Speciesn{11+q2});RO11=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO11=horzcat('kine',num2str(n));Rkonst11=horzcat('K',num2str(n));
       eval([sprintf(RO11) '=addreaction(SMV2,sprintf(R11));']); 
       eval([sprintf(KLO11) '= addkineticlaw(eval(sprintf(RO11)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst11),Parameters.dimparams.MonAssNuc(q2),V,U2);']); 
       set(eval(sprintf(KLO11)),'ParameterVariableName', {sprintf(Rkonst11)});n=n+1;
       %Dissociation  of Dimers in Cytoplasm
       R12=horzcat(Species{11+q2},' ','->',' ',Species{11+q2}(1:4),' ','+',' ',Species{11+q2}(5:end));RO12=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO12=horzcat('kine',num2str(n));Rkonst12=horzcat('K',num2str(n));
       eval([sprintf(RO12) '=addreaction(SMV2,sprintf(R12));']); 
       eval([sprintf(KLO12) '= addkineticlaw(eval(sprintf(RO12)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst12),Parameters.dimparams.DimDissCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO12)),'ParameterVariableName', {sprintf(Rkonst12)});n=n+1;
       %Dissociation  of Dimers in Nucleus
       R13=horzcat(Speciesn{11+q2},' ','->',' ',Speciesn{11+q2}(1:4),'n',' ','+',' ',Speciesn{11+q2}(5:end));RO13=horzcat('Nuc_Ex_Mon',Species{11+q2},num2str(n));KLO13=horzcat('kine',num2str(n));Rkonst13=horzcat('K',num2str(n));
       eval([sprintf(RO13) '=addreaction(SMV2,sprintf(R13));']); 
       eval([sprintf(KLO13) '= addkineticlaw(eval(sprintf(RO13)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst13),Parameters.dimparams.DimDissNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO13)),'ParameterVariableName', {sprintf(Rkonst13)});n=n+1;
       %Nuclear Import
       R14=horzcat(Species{11+q2},' ','->',' ',Speciesn{11+q2});RO14=horzcat('Nuc_Im_Dim',Species{q2+11},num2str(n));KLO14=horzcat('kine',num2str(n));Rkonst14=horzcat('K',num2str(n));
       eval([sprintf(RO14) '=addreaction(SMV2,sprintf(R14));']); 
       eval([sprintf(KLO14) '= addkineticlaw(eval(sprintf(RO14)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst14),Parameters.dimparams.NucImp(q2),V,U1);']); 
       set(eval(sprintf(KLO14)),'ParameterVariableName', {sprintf(Rkonst14)});n=n+1;
       %Nuclear Export
       R15=horzcat(Speciesn{11+q2},' ','->',' ',Species{11+q2});RO15=horzcat('Nuc_Ex_Dim',Species{q2+11},num2str(n));KLO15=horzcat('kine',num2str(n));Rkonst15=horzcat('K',num2str(n));
       eval([sprintf(RO15) '=addreaction(SMV2,sprintf(R15));']); 
       eval([sprintf(KLO15) '= addkineticlaw(eval(sprintf(RO15)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst15),Parameters.dimparams.NucExp(q2),V,U1);']); 
       set(eval(sprintf(KLO15)),'ParameterVariableName', {sprintf(Rkonst15)});n=n+1;
       %Free Dimer Degradation in Cytoplasm
       R16=horzcat(Species{11+q2},' ','->',' ','null');RO16=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO16=horzcat('kine',num2str(n));Rkonst16=horzcat('K',num2str(n));
       eval([sprintf(RO16) '=addreaction(SMV2,sprintf(R16));']); 
       eval([sprintf(KLO16) '= addkineticlaw(eval(sprintf(RO16)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst16),Parameters.dimparams.DimDegCyt(q2),V,U1);']); 
       set(eval(sprintf(KLO16)),'ParameterVariableName', {sprintf(Rkonst16)});n=n+1;
       %Free Dimer Degradation in Nucleus
       R17=horzcat(Speciesn{11+q2},' ','->',' ','null');RO17=horzcat('Cyto_Deg_Dim',Species{q2+11},num2str(n));KLO17=horzcat('kine',num2str(n));Rkonst17=horzcat('K',num2str(n));
       eval([sprintf(RO17) '=addreaction(SMV2,sprintf(R17));']); 
       eval([sprintf(KLO17) '= addkineticlaw(eval(sprintf(RO17)),sprintf(M));']);
       eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst17),Parameters.dimparams.DimDegNuc(q2),V,U1);']); 
       set(eval(sprintf(KLO17)),'ParameterVariableName', {sprintf(Rkonst17)});n=n+1;
end
%========================================================
                                                     
for q3=1:16
        %Association of trimers in Cytoplasm
        R18=horzcat(Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end),' ','->',' ',Species{q3+19});RO18=horzcat('Trim_Ass_Cyto',Species{q3+19},num2str(n));KLO18=horzcat('kine',num2str(n));Rkonst18=horzcat('K',num2str(n));
        eval([sprintf(RO18) '=addreaction(SMV2,sprintf(R18));']); 
        eval([sprintf(KLO18) '= addkineticlaw(eval(sprintf(RO18)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst18),Parameters.trimparams.Trim_Ass_Cyt(q3),V,U2);']); 
        set(eval(sprintf(KLO18)),'ParameterVariableName', {sprintf(Rkonst18)});n=n+1;
        %Association of trimers in Nucleus
        R19=horzcat(Speciesn{q3+19}(1:4),'n',' ','+',' ',Speciesn{q3+19}(5:end),' ','->',' ',Speciesn{q3+19});RO19=horzcat('Trim_Ass_Nuc',Species{q3+19},num2str(n));KLO19=horzcat('kine',num2str(n));Rkonst19=horzcat('K',num2str(n));
        eval([sprintf(RO19) '=addreaction(SMV2,sprintf(R19));']); 
        eval([sprintf(KLO19) '= addkineticlaw(eval(sprintf(RO19)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst19),Parameters.trimparams.Trim_Ass_Nuc(q3),V,U2);']); 
        set(eval(sprintf(KLO19)),'ParameterVariableName', {sprintf(Rkonst19)});n=n+1;
        %Dissociation of trimers in Cytoplasm
        R20=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4),' ','+',' ',Species{q3+19}(5:end));RO20=horzcat('Trim_Dis_Cyt',Species{q3+19},num2str(n));KLO20=horzcat('kine',num2str(n));Rkonst20=horzcat('K',num2str(n));
        eval([sprintf(RO20) '=addreaction(SMV2,sprintf(R20));']); 
        eval([sprintf(KLO20) '= addkineticlaw(eval(sprintf(RO20)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst20),Parameters.trimparams.Trim_Dis_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO20)),'ParameterVariableName', {sprintf(Rkonst20)});n=n+1;
        %Dissociation of trimers in Nucleus
        R21=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n',' ','+',' ',Species{q3+19}(5:end),'n');RO21=horzcat('Trim_Dis_Nuc',Species{q3+19},num2str(n));KLO21=horzcat('kine',num2str(n));Rkonst21=horzcat('K',num2str(n));
        eval([sprintf(RO21) '=addreaction(SMV2,sprintf(R21));']); 
        eval([sprintf(KLO21) '= addkineticlaw(eval(sprintf(RO21)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst21),Parameters.trimparams.Trim_Dis_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO21)),'ParameterVariableName', {sprintf(Rkonst21)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Cytoplasm
        R22=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(5:end));RO22=horzcat('Const_Deg_Trim_Ikb_Cyt',Species{q3+19},num2str(n));KLO22=horzcat('kine',num2str(n));Rkonst22=horzcat('K',num2str(n));
        eval([sprintf(RO22) '=addreaction(SMV2,sprintf(R22));']); 
        eval([sprintf(KLO22) '= addkineticlaw(eval(sprintf(RO22)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst22),Parameters.trimparams.Const_Trim_Deg_IkB_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO22)),'ParameterVariableName', {sprintf(Rkonst22)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> NFkB in Nucleus
        R23=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(5:end),'n');RO23=horzcat('Const_Deg_Trim_Ikb_Nuc',Species{q3+19},num2str(n));KLO23=horzcat('kine',num2str(n));Rkonst23=horzcat('K',num2str(n));
        eval([sprintf(RO23) '=addreaction(SMV2,sprintf(R23));']); 
        eval([sprintf(KLO23) '= addkineticlaw(eval(sprintf(RO23)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst23),Parameters.trimparams.Const_Trim_Deg_IkB_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO23)),'ParameterVariableName', {sprintf(Rkonst23)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Cytoplasm
        R24=horzcat(Species{q3+19},' ','->',' ',Species{q3+19}(1:4));RO24=horzcat('Const_Deg_Trim_NFkb_Cyt',Species{q3+19},num2str(n));KLO24=horzcat('kine',num2str(n));Rkonst24=horzcat('K',num2str(n));
        eval([sprintf(RO24) '=addreaction(SMV2,sprintf(R24));']); 
        eval([sprintf(KLO24) '= addkineticlaw(eval(sprintf(RO24)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst24),Parameters.trimparams.Const_Trim_Deg_NFkb_Cyt(q3),V,U1);']); 
        set(eval(sprintf(KLO24)),'ParameterVariableName', {sprintf(Rkonst24)});n=n+1;
        %Constitutive Degradation of trimer : Ikb:NFkB -> Ikb in Nucleus
        R25=horzcat(Speciesn{q3+19},' ','->',' ',Species{q3+19}(1:4),'n');RO25=horzcat('Const_Deg_Trim_NFkb_Nuc',Species{q3+19},num2str(n));KLO25=horzcat('kine',num2str(n));Rkonst25=horzcat('K',num2str(n));
        eval([sprintf(RO25) '=addreaction(SMV2,sprintf(R25));']); 
        eval([sprintf(KLO25) '= addkineticlaw(eval(sprintf(RO25)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst25),Parameters.trimparams.Const_Trim_Deg_NFkb_Nuc(q3),V,U1);']); 
        set(eval(sprintf(KLO25)),'ParameterVariableName', {sprintf(Rkonst25)});n=n+1;
        %NEMO Mediated Degradation of Ikb in Trimer in Cytoplasm
        R26=horzcat(Species{q3+19},' ','+',' ','NEMO_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO26=horzcat('NEMO_DEg_inTrim',Species{q3+19},num2str(n));KLO26=horzcat('kine',num2str(n));Rkonst26=horzcat('K',num2str(n));
        eval([sprintf(RO26) '=addreaction(SMV2,sprintf(R26));']); 
        eval([sprintf(KLO26) '= addkineticlaw(eval(sprintf(RO26)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst26),Parameters.trimparams.NEMO_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO26)),'ParameterVariableName', {sprintf(Rkonst26)});n=n+1;
        %NIK Mediated Degradation of Ikb in Trimer in Cytoplasm
        R27=horzcat(Species{q3+19},' ','+',' ','NIK_ACTIVITY',' ','->',' ',Species{q3+19}(5:end));RO27=horzcat('NIK_DEg_inTrim',Species{q3+19},num2str(n));KLO27=horzcat('kine',num2str(n));Rkonst27=horzcat('K',num2str(n));
        eval([sprintf(RO27) '=addreaction(SMV2,sprintf(R27));']); 
        eval([sprintf(KLO27) '= addkineticlaw(eval(sprintf(RO27)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst27),Parameters.trimparams.NIK_Deg_Ikb_inTrim(q3),V,U2);']); 
        set(eval(sprintf(KLO27)),'ParameterVariableName', {sprintf(Rkonst27)});n=n+1;
        %Nuclear Import of Trimers
        R28=horzcat(Species{q3+18},' ','->',Speciesn{q3+18});RO28=horzcat('Nuc_Imp_Trim',Species{q3+18},num2str(n));KLO28=horzcat('kine',num2str(n));Rkonst28=horzcat('K',num2str(n));
        eval([sprintf(RO28) '=addreaction(SMV2,sprintf(R28));']); 
        eval([sprintf(KLO28) '= addkineticlaw(eval(sprintf(RO28)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst28),Parameters.trimparams.Nuc_Imp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO28)),'ParameterVariableName', {sprintf(Rkonst28)});n=n+1;
        %Nuclear Export of Trimers
        R29=horzcat(Speciesn{q3+19},' ','->',Species{q3+19});RO29=horzcat('Nuc_Exp_Trip',Species{q3+19},num2str(n));KLO29=horzcat('kine',num2str(n));Rkonst29=horzcat('K',num2str(n));
        eval([sprintf(RO29) '=addreaction(SMV2,sprintf(R29));']); 
        eval([sprintf(KLO29) '= addkineticlaw(eval(sprintf(RO29)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst29),Parameters.trimparams.Nuc_Exp_Trim(q3),V,U1);']); 
        set(eval(sprintf(KLO29)),'ParameterVariableName', {sprintf(Rkonst29)});n=n+1;
end
%=========================================================
%Special Reactions
   Parameters.p100_Ass            = 9e-2;     % p100 association rate to form IkBd in cytoplasm %%% 1.2E-3
   %                           K253 
   Parameters.Ikbd_Diss           = 1.2e-5;   % IkBd dissocation rate to form p100 in cytoplasm %%% 1.2E-2
   %                           K254
   Parameters.NIK_mediated_p52    = 3.5e-3;   % NIK mediated generation of P52 through processing
   %                           K255
   
        R30=horzcat('2',' ',Species{5},' ','->',Species{4});RO30=horzcat('p100_Ass',Species{5},num2str(n));KLO30=horzcat('kine',num2str(n));Rkonst30=horzcat('K',num2str(n));
        eval([sprintf(RO30) '=addreaction(SMV2,sprintf(R30));']); 
        eval([sprintf(KLO30) '= addkineticlaw(eval(sprintf(RO30)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst30),Parameters.p100_Ass,V,U2);']); 
        set(eval(sprintf(KLO30)),'ParameterVariableName', {sprintf(Rkonst30)});n=n+1;

        R31=horzcat(Species{4},' ','->',' ','2',' ',Species{5});RO31=horzcat('p100_Ass',Species{4},num2str(n));KLO31=horzcat('kine',num2str(n));Rkonst31=horzcat('K',num2str(n));
        eval([sprintf(RO31) '=addreaction(SMV2,sprintf(R31));']); 
        eval([sprintf(KLO31) '= addkineticlaw(eval(sprintf(RO31)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst31),Parameters.Ikbd_Diss,V,U1);']); 
        set(eval(sprintf(KLO31)),'ParameterVariableName', {sprintf(Rkonst31)});n=n+1;

        R32=horzcat(Species{5},' ','+',' ',Species{end},' ','->',' ',Species{10});RO32=horzcat('p52_processing',Species{5},num2str(n));KLO32=horzcat('kine',num2str(n));Rkonst32=horzcat('K',num2str(n));
        eval([sprintf(RO32) '=addreaction(SMV2,sprintf(R32));']); 
        eval([sprintf(KLO32) '= addkineticlaw(eval(sprintf(RO32)),sprintf(M));']);
        eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rkonst32),Parameters.NIK_mediated_p52,V,U2);']); 
        set(eval(sprintf(KLO32)),'ParameterVariableName', {sprintf(Rkonst32)});n=n+1;
%===========================================================
%Chequer Boxes of Inducibility 
%=============================
%                   
%Add Rules
%==========

RfntIkba = 'tIkba = ( Kct1 * ( ( 1 + ( W11 * ( RelAp50n / Kd11 ) ^ 3 ) + ( W12 * ( RelAp52n / Kd12 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd11 ) ^ 3 + ( RelAp52n / Kd12 ) ^ 3 ) ) ) - Kmdeg1 * tIkba';
              
RfntIkbb = 'tIkbb = ( Kct2 * ( ( 1 + ( W21 * ( RelAp50n / Kd21 ) ^ 3 ) + ( W22 * ( RelAp52n / Kd22 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd21 ) ^ 3 + ( RelAp52n / Kd22 ) ^ 3 ) ) ^ Kdel1 ) - Kmdeg3 * tIkbb';
                 
RfntIkbe = 'tIkbe = ( Kct3 * ( ( 1 + ( W31 * ( RelAp50n / Kd31 ) ^ 3 ) + ( W32 * ( RelAp52n / Kd32 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd31 ) ^ 3 + ( RelAp52n / Kd32 ) ^ 3 ) ) ^ Kdel2 ) - Kmdeg3 * tIkbe';
                
Rfntp100 = 'tp100 = ( Kct5 * ( ( 1 + ( W41 * ( RelAp50n / Kd41 ) ^ 3 ) + ( W42 * ( RelAp52n / Kd42 ) ^ 3 ) ) / ( 1 + ( RelAp50n / Kd41 ) ^ 3 + ( RelAp52n / Kd42 ) ^ 3 ) ) ^ Kdel3 ) - Kmdeg5 * tp100';
                 
RfntRelA = 'tRelA = Kct7  - Kmdeg7 * tRelA';

RfntRelB = 'tRelB =  Kct8 - Kmdeg8 * tRelB';

RfntR_lB = 'tR_lB =  Kct11  - Kmdeg11 * tR_lB';
                   
Rfntp50  = 'tp50  = Kct9 - Kmdeg9 * tp50';
for krow= 1:size(Parameters.Kd,1)
  for kcolumn= 1:size(Parameters.Kd,2)
      Rconst=horzcat('Kd',num2str(krow),num2str(kcolumn)); Rconst1=horzcat('W',num2str(krow),num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),Parameters.Kd(krow,kcolumn),V,U);']);n=n+1;
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst1),Parameters.W(krow,kcolumn));']);n=n+1; 
  end
end
%============================================
  RlO1 = addrule(SMV2,sprintf(RfntIkba),'rate');
  RlO2 = addrule(SMV2,sprintf(RfntIkbb),'rate');
  RlO3 = addrule(SMV2,sprintf(RfntIkbe),'rate');
  RlO4 = addrule(SMV2,sprintf(Rfntp100),'rate');
  RlO5 = addrule(SMV2,sprintf(RfntRelA),'rate');
  RlO6 = addrule(SMV2,sprintf(RfntRelB),'rate');
  RlO7 = addrule(SMV2,sprintf(RfntR_lB),'rate');
  RlO8 = addrule(SMV2,sprintf(Rfntp50),'rate');
%============================================
Ca='ConstantValue';%VA='Value'

%DELAY
%=====
%kdels are 1 during euilibration
for kcolumn= 1:5
      Rconst=horzcat('Kdel',num2str(kcolumn));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(Rconst),1,Ca,false);']);
      n=n+1; 
end
%kdels are 0 at stimulation 0 min.
for P=1:4
 E=horzcat('Eee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(0));
                                         T=horzcat('time >= ',num2str(10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%kdels are 1 at stimulation phase once delay is over.
DELAY=[37 37 45 0]; %beta epsilon p100 RelB 
for P=1:4
 E=horzcat('Eeee',num2str(P));N=horzcat('Kdel',num2str(P),' ','=',' ',num2str(1));
                                         T=horzcat('time >= ',num2str(DELAY(P)+10000));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end
%RelB star issue
%===============

      RconsT=horzcat('Kct',num2str(8));RconsTd=horzcat('Kmdeg',num2str(8));
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsT),Parameters.monparams.Const_trnx(8),Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Param',num2str(n))) '= addparameter(SMV2,sprintf(RconsTd),Parameters.monparams.mRNA_Deg(8),Ca,false);']); n=n+1; 
      EEE1=addevent(SMV2,'time >= 10000','Kct8 = 0');%EEE2=addevent(SMV2,'time >= 150000','Kmdeg8 = 0');

      RconsTT=horzcat('Kct',num2str(11));RconsTTd=horzcat('Kmdeg',num2str(11));
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTT),0,Ca,false);']); n=n+1; 
      eval([sprintf(horzcat('Paramt',num2str(n))) '= addparameter(SMV2,sprintf(RconsTTd),0,Ca,false);']); n=n+1; 
      %EEE3=addevent(SMV2,'time >= 150000','Kct11 =  3.6e-5');EEE4=addevent(SMV2,'time >= 150000','Kmdeg11 = 4.6e-3');
      E=horzcat('EEE',num2str(3));N=horzcat('Kct',num2str(11),' ','=',' ',num2str(Parameters.monparams.Const_trnx(11))); T=horzcat('time >= ',num2str(10000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);       
      E=horzcat('EEE',num2str(4));N=horzcat('Kmdeg',num2str(11),' ','=',' ',num2str(Parameters.monparams.mRNA_Deg(11))); %T=horzcat('time >= ',num2str(150000));
      eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);                                                                                                  
%Effective RelA heterodimer concentration in nucleus involved with RelB synthesis
evalues =  [1 .01];etime = [0 60];
 for Q=0:1:60
 Kdel5 = interp1(etime, evalues, Q,'pchip');E=horzcat('Eeff',num2str(Q));N=horzcat('Kdel5',' ','=',' ',num2str(Kdel5));
                                         T=horzcat('time >= ',num2str(10030+Q));%RelA effective concentration start going down after 30 mins
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end 
%SIGNAL
%======
%TNFc
%cvalues = [1 60 100 65 50 36 20 10 3];ctime = [0 5 10 15 20 25 30 360 480];
%cvalues =  [1 60 100 65 50 36 20 10 3];ctime = [0 5 10 15 20 25 30 360 480];
%======================================
%Add Parameters of Inducibility
%==============================
for Q=0:1:480
 NEMO = interp1(ptime, pvalues, Q,'pchip');E=horzcat('E',num2str(Q));N=horzcat('NEMO_ACTIVITY',' ','=',' ',num2str(NEMO));
                                         T=horzcat('time >= ',num2str(10000+Q));
 eval([sprintf(E) '=addevent(SMV2,sprintf(T),sprintf(N));']);
end

%SIMULATE MODEL
%==============
CnfgstO = addconfigset(SMV2, 'De3');
set(CnfgstO.CompileOptions, 'UnitConversion', false);
set(CnfgstO.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO,'TimeUnits',' ');
set(CnfgstO,'SolverType','sundials');
set(CnfgstO,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO)
[t11,x11,names11]= sbiosimulate(SMV2,CnfgstO);
%ADD VARIANT 1
%==============
estvarObj1 = addvariant (SMV2, 'nfkb2kn');
addcontent(estvarObj1, {'parameter','K34', 'Value',0});
CnfgstO1 = addconfigset(SMV2, 'nfkb2kn');
set(CnfgstO1.CompileOptions, 'UnitConversion', false);
set(CnfgstO1.CompileOptions, 'DimensionalAnalysis', false);
%set(CnfgstO.CompileOptions, 'DefaultSpeciesDimension', 'substance');
%set(CnfgstO1,'TimeUnits',' ');
set(CnfgstO1,'SolverType','sundials');
set(CnfgstO1,'StopTime', 10480);
sbioaccelerate(SMV2,CnfgstO1,estvarObj1)
[t12,x12,names12]= sbiosimulate(SMV2,CnfgstO1,estvarObj1);

%PLOT RESULTS
%============
figure('NumberTitle','on','units','normalized','outerposition',[0 0 1 1])
f=.013;Time1=10000;Time2=10480;

for a=0:1:11
 
 tt=horzcat('t',num2str(a+1));
 xx=horzcat('x',num2str(a+1));
 T=eval(tt);
 X=eval(xx);

%Free Absolute Nuclear NFkB
if a <= 5
subplot_tight(6,4,a*4+1,f)
end
if a >= 6
subplot_tight(6,4,(a-6)*4+3,f)
end
plot(T,X(:,24)+X(:,26),'r','linewidth',2);hold on   %RelAp50n
plot(T,X(:,28)+X(:,36),'g','linewidth',3);hold on  %RelBp50n+R_lBp50n   
plot(T,X(:,28),'c','linewidth',1);hold on %RelBp50n 
plot(T,X(:,36),'b','linewidth',1);hold off %R_lBp50n
if a==0
legend('RelAp50n','RelBp50nT','RelBp50n','R_lBp50n')
legend('Location','NorthEast')
end
xlim([Time1 Time2])
ylim([0 150])
set(gca,'xtick',[0 30 60 90 120 150 180 210 240 270 300 330 360 390 420 450 480]+10000)
set(gca,'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' '])
set(gca,'XGrid','on')

%Transcript fold induction
for D=1:size(T)
    if T(D)>= 10000
       break
    end
end
if a <= 5
subplot_tight(6,4,a*4+2,f)
end
if a >= 6
subplot_tight(6,4,(a-6)*4+4,f)
end
plot(T,X(:,71)/X(D,71),'b');hold on   %tIkba
plot(T,X(:,72)/X(D,72),'c');hold on   %tIkbb
plot(T,X(:,73)/X(D,73),'m');hold on   %tIkbe
plot(T,X(:,75)/X(D,75),'k');hold on   %tp100
plot(T,X(:,77)/X(D,77),'r');hold on   %tRelA
plot(T,((X(:,78)+X(:,81))/(X(D,78)+X(D,81))),'g');hold on   %tRelB
plot(T,X(:,79)/X(D,79),'k--');hold on  %tp50
if a==0
legend('tIkBa','tIkBb','tIkBe','tp100','tRElA','tRelB','tp50')
legend('Location','NorthEast')
end
xlim([Time1 Time2])
ylim([0 15])
set(gca,'xtick',[0 30 60 90 120 150 180 210 240 270 300 330 360 390 420 450 480]+10000)
set(gca,'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' '])
set(gca,'XGrid','on')

end


clearvars -except H t1 x1 names1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 SMV2
elapsedtime=toc;H=H+1;
fprintf('Simulation no. %-5.2f : elapsed time : %-5.2f', H, elapsedtime);
