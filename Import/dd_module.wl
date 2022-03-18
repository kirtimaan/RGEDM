(* ::Package:: *)

(*Image Cropping Module*)
ImCrop[image_]:=Module[{nimage},
nimage=Rasterize[image,"Image",ImageResolution->300]];
(*Export["figname.png",ImCrop[pl1]]*)


(* ::Subsubsection::Closed:: *)
(*Numerical Input and Wilson Coeffs*)


(*Numerical Input*)
(*For light quarks we sue current quakr masses from PDG ealuated in the MSbar scheme at \[Mu]=2GeV*)
mt=173.0;
mb=4.2;
mc=1.3;
ms=95.0/1000.0;
md=4.70/1000.0;
mu=2.2/1000.0;
mZ=91.188;
\[Alpha]smZ=\!\(TraditionalForm\`0.11848964091136216\);
gsmZ=Sqrt[\[Alpha]smZ*4*Pi];
mp=0.938272046;
(*Xenon mass number and atomic number*)
Axe=131; 
Zxe = 54;
(**Nuclear parameters**)
(*Defined at mu = 2 GeV*)
fTq=0.016;
fTG=0.80;(*NNLO*)
Q2=0.346;
(*(*CT14NNLO pdf set Q=1.3 GeV*)
G2=0.388;*)
(*CT14NNLO pdf set Q=2 GeV*)
G2=0.4159;
(*{u,d,s}*)
(*(*CT14NNLO pdf set Q=1.3 GeV*)
Q2list={0.3796,0.2029,0.0295,0.0,0.0,0.0};*)
(*CT14NNLO pdf set Q=2 GeV*)
Q2list={0.3481,0.1902,0.0352,0.0107,0.0,0.0};
(*CTEQ pdf evaluated at Z mass*)
Q2listZ={0.254,0.146,0.052,0.038,0.024,0.0};
(*Define for neutron assuming isospin sym*)

Q2listn={0.1902,0.3481,0.0352,0.0107,0.0,0.0};
(*CTEQ pdf evaluated at Z mass*)
Q2listnZ={0.146,0.254,0.052,0.038,0.024,0.0};
fTqlist={0.018,0.030,0.043,0.0,0.0,0.0};
fTqlistn={0.015,0.034,0.043,0.0,0.0,0.0};
\[CapitalDelta]up=0.84;
\[CapitalDelta]dp=-0.43;
\[CapitalDelta]sp=-0.09;
\[CapitalDelta]un=\[CapitalDelta]dp;
\[CapitalDelta]dn=\[CapitalDelta]up;
\[CapitalDelta]sn=\[CapitalDelta]sp;

mmL={mu,md,ms,mc,mb,mt};

(*Scale upto which we run*)
mul=2.0;



\[Alpha]srun2[\[Mu]i_,\[Mu]f_,nf_,L_,gsL_]:=Module[{\[Beta]0=11 -2/3nf,\[Beta]1= 102 - 38/3 nf , \[Beta]2= 2857/2 -5033/18 nf + 325/54 nf^2,\[Beta]3= 29243-6946.30 nf + 405.089 nf^2 + 1.49931 nf^3,ret=0},
ret=NDSolveValue[{g'[x]== g[x](-\[Beta]0 (g[x]^2/16/Pi^2)  -\[Beta]1 (g[x]^2/16/Pi^2)^2   -\[Beta]2 (g[x]^2/16/Pi^2)^3   -\[Beta]3 (g[x]^2/16/Pi^2)^4),g[Log[L]]==gsL},g,{x,Log[\[Mu]i],Log[\[Mu]f]}] ]

(*t1=\[Alpha]srun2[mt,mb,5,mZ,gsmZ];*)
t1=\[Alpha]srun2[mt,mb,5,mZ,\!\(TraditionalForm\`gsmZ\)];
t0=\[Alpha]srun2[10000,mt,6,mt,t1[Log[mt]]];
t2=\[Alpha]srun2[mb,mc,4,mb,t1[Log[mb]]];
t3=\[Alpha]srun2[mc,mp,3,mc,t2[Log[mc]]];
(*\[Alpha]srun[x_]:=Module[{},Which[x>=mt,t0[Log[x]],x<mt&&x\[GreaterEqual]mb,t1[Log[x]],x<mb &&x\[GreaterEqual] mc,t2[Log[x]],x<mc && x\[GreaterEqual] mp,t3[Log[x]],x<mp ,t3[Log[mp]]]];*)
\[Alpha]srun[x_]:=Module[{},Which[x>=mt,t0[Log[x]]^2/(4 Pi),x<mt&&x>=mb,t1[Log[x]]^2/(4 Pi),x<mb &&x>= mc,t2[Log[x]]^2/(4 Pi),x<mc && x>= mp,t3[Log[x]]^2/(4 Pi),x<mp ,t3[Log[mp]]^2/(4 Pi)]];



(*Definition of Evolution Matrix*)

\[CapitalLambda]nf[nf_]:=Module[{ret=0},Switch[nf,5,ret=0.213066,4,ret=0.297608,3,ret=0.339872,_?NumericQ,ret=0,_,ret= \[CapitalLambda][nf]];ret];

\[Alpha]srunapprox[\[Mu]_,nf_]:=Module[{\[Beta]0=11 -2/3nf,\[Beta]1= 102 - 38/3 nf , \[Beta]2= 2857/2 -5033/18 nf + 325/54 nf^2,\[Beta]3= 29243-6946.30 nf + 405.089 nf^2 + 1.49931 nf^3,t=Log[\[Mu]^2/(\[CapitalLambda]nf[nf]^2)],ret=0},
ret=(4 Pi/(\[Beta]0 t))*(
1 - \[Beta]1* Log[t]/(\[Beta]0^2 t)
+ \[Beta]1^2/(\[Beta]0^4 t^2)* ( (Log[t] -1/2)^2 - 5/4 + \[Beta]0 \[Beta]2/\[Beta]1^2)
-1/(\[Beta]0^6 t^3)*(\[Beta]1^3 (Log[t]^3 -5/2 Log[t]^2 - 2 Log[t] + 1/2) + 3 \[Beta]0 \[Beta]1 \[Beta]2 - 1/2 \[Beta]0^2  \[Beta]3 )
)  ];

\[Beta]tilde[\[Mu]_,nf_]:=Module[{\[Beta]0=11 -2/3nf,\[Beta]1= 102 - 38/3 nf , \[Beta]2= 2857/2 -5033/18 nf + 325/54 nf^2,\[Beta]3= 29243-6946.30 nf + 405.089 nf^2 + 1.49931 nf^3,ret=0,x=\[Alpha]srun[\[Mu]]/(4 Pi)},ret=-\[Beta]0*x - \[Beta]1*x^2 - \[Beta]2 x^3 -\[Beta]3 x^4];

\[Gamma]mrun[\[Mu]_,nf_]:=Module[{\[Gamma]0=8,\[Gamma]1=404/3 - 40/9 nf,\[Gamma]2=2498 - (4432/27 + 320/3 Zeta[3])nf - 280/81 nf^2,\[Gamma]3=50659 - 9783.04 nf + 141.395 nf^2 + 2.96613 nf^3 ,x=\[Alpha]srun[\[Mu]]/(4 Pi),ret=0 },ret= -\[Gamma]0 * x -\[Gamma]1 * x^2 - \[Gamma]2*x^3 -\[Gamma]3*x^4]


r[t_,nf_,\[Mu]l_,\[Mu]h_]:=Module[{\[Beta]0= 11 -2/3 nf},(\[Alpha]srun[\[Mu]l]/\[Alpha]srun[\[Mu]h])^(-(2/(3 \[Beta]0)) (16/3 +  t))];

R0qgScalar[nf_,\[Mu]l_,\[Mu]h_]:=Module[{Rqq,Rqqp,Rqg,Rgq,Rgg,JJ= ConstantArray[1,{nf,nf}],qqmat,Qgmat,gQmat,ret},Rqq=1;Rqqp=0;Rqg=2*(\[Gamma]mrun[\[Mu]h,nf]- \[Gamma]mrun[\[Mu]l,nf])/\[Beta]tilde[\[Mu]h,nf]; Rgq=0;Rgg= \[Beta]tilde[\[Mu]l,nf]/\[Beta]tilde[\[Mu]h,nf]; qqmat=(Rqq -Rqqp)IdentityMatrix[nf] + Rqqp*JJ; gQmat= ConstantArray[Rgq,{nf}];Qgmat=Join[ConstantArray[Rqg,{nf,1}],{{Rgg}}];ret=  Join[Join[qqmat,{gQmat}], Qgmat,2]];

R0qgSpin2[nf_,\[Mu]l_,\[Mu]h_]:=Module[{RqqI,Rqqp,Rqg,Rgq,Rgg,JJ= ConstantArray[1,{nf,nf}],qqmat,Qgmat,gQmat,ret,rr= r[nf,nf,\[Mu]l,\[Mu]h],rr0= r[0,nf,\[Mu]l,\[Mu]h]},RqqI=rr0;Rqqp=1/nf ((16 rr + 3 nf )/(16 + 3 nf) - rr0);Rqg=16(1 -rr)/(16 + 3 nf); Rgq=3*(1 -rr)/(16 + 3 nf);Rgg= (16 + 3 nf * rr)/(16 + 3 nf); qqmat=RqqI*IdentityMatrix[nf] + Rqqp*JJ; gQmat= ConstantArray[Rgq,{nf}];Qgmat=Join[ConstantArray[Rqg,{nf,1}],{{Rgg}}];ret=  Join[Join[qqmat,{gQmat}], Qgmat,2]];



(*Definition of matching matrix*)
M0qgScalar[nf_,\[Mu]Q_,mQ_]:=Module[{MgQ,Mgg,qqmat,Qgmat,gQmat,ret, \[Alpha]sp= \[Alpha]srun[\[Mu]Q]},Mgg=1 - (\[Alpha]sp/(3 Pi))Log[\[Mu]Q/mQ];MgQ=-(\[Alpha]sp/(12 Pi))(1 +\[Alpha]sp/(4 Pi) *(11 - 4/3 Log[\[Mu]Q/mQ]) ); qqmat=Join[IdentityMatrix[nf-1],ConstantArray[{0,0},nf-1],2] ; gQmat= Join[ConstantArray[0,{nf-1}],{MgQ,Mgg}];ret=  Join[qqmat, {gQmat}]];

M0qgSpin2[nf_,\[Mu]Q_,mQ_]:=Module[{MgQ,Mgg,qqmat,Qgmat,gQmat,ret, \[Alpha]sp= \[Alpha]srun[\[Mu]Q]},Mgg=1 ;MgQ=-(\[Alpha]sp/(3 Pi))Log[\[Mu]Q/mQ]; qqmat=Join[IdentityMatrix[nf-1],ConstantArray[{0,0},nf-1],2] ; gQmat= Join[ConstantArray[0,{nf-1}],{MgQ,Mgg}];ret=  Join[qqmat, {gQmat}]];



(*RGE of wilson coefficient*)
cgq0={cu0,cd0,cs0,cc0,cb0,ct0,cg0};
cgq2={cu2,cd2,cs2,cc2,cb2,ct0,cg2};
renc0[cgq0_,\[Mu]0_]:=Module[{ret},ret=R0qgScalar[3,\[Mu]0,mc] . M0qgScalar[4,mc,mc] . R0qgScalar[4,mc,mb] . M0qgScalar[5,mb,mb] . R0qgScalar[5,mb,mt] . M0qgScalar[6,mt,mt] . R0qgScalar[6,mt,mt] . cgq0 ];
renc2[cgq2_,\[Mu]0_]:=Module[{ret},ret=R0qgSpin2[3,\[Mu]0,mc] . M0qgSpin2[4,mc,mc] . R0qgSpin2[4,mc,mb] . M0qgSpin2[5,mb,mb] . R0qgSpin2[5,mb,mt] . M0qgSpin2[6,mt,mt] . R0qgSpin2[6,mt,mt] . cgq2 ];

(*Define with arbitrary mu*)
renc0mu[cgq2_,\[Mu]0_]:=Catch[Module[{ret},If[Length[cgq2]!=7,Message["Bad argument given to function renc2mu: ",cgq2];Throw[$Failed]];Which[\[Mu]0>=mt,ret=cgq2,\[Mu]0< mt && \[Mu]0>=  mb,ret=R0qgScalar[5,\[Mu]0,mt] . M0qgScalar[6,mt,mt] . R0qgScalar[6,mt,mt] . cgq2,\[Mu]0< mb && \[Mu]0>=  mc,ret=R0qgScalar[4,mc,mb] . M0qgScalar[5,mb,mb] . R0qgScalar[5,mb,mt] . M0qgScalar[6,mt,mt] . R0qgScalar[6,mt,mt] . cgq2,\[Mu]0< mc && \[Mu]0>=  mp,ret=R0qgScalar[3,\[Mu]0,mc] . M0qgScalar[4,mc,mc] . R0qgScalar[4,mc,mb] . M0qgScalar[5,mb,mb] . R0qgScalar[5,mb,mt] . M0qgScalar[6,mt,mt] . R0qgScalar[6,mt,mt] . cgq2,\[Mu]0<mp,Throw[$Failed,\[Mu]0]]]];
renc2mu[cgq2_,\[Mu]0_]:=Catch[Module[{ret},If[Length[cgq2]!=7,Message["Bad argument given to function renc2mu: ",cgq2];Throw[$Failed]];Which[\[Mu]0>=mt,ret=cgq2,\[Mu]0< mt && \[Mu]0>=  mb,ret=R0qgSpin2[5,\[Mu]0,mt] . M0qgSpin2[6,mt,mt] . R0qgSpin2[6,mt,mt] . cgq2,\[Mu]0< mb && \[Mu]0>=  mc,ret=R0qgSpin2[4,mc,mb] . M0qgSpin2[5,mb,mb] . R0qgSpin2[5,mb,mt] . M0qgSpin2[6,mt,mt] . R0qgSpin2[6,mt,mt] . cgq2,\[Mu]0< mc && \[Mu]0>=  mp,ret=R0qgSpin2[3,\[Mu]0,mc] . M0qgSpin2[4,mc,mc] . R0qgSpin2[4,mc,mb] . M0qgSpin2[5,mb,mb] . R0qgSpin2[5,mb,mt] . M0qgSpin2[6,mt,mt] . R0qgSpin2[6,mt,mt] . cgq2,\[Mu]0<mp,Throw[$Failed,\[Mu]0]]]];



(*Wilson Coeffs*)
fG[M_,m_,mchi_,Alfas_,gDM_]:=Alfas*gDM^2*(mchi*(-m^10-(M^2-mchi^2)^4*(2*M^2-mchi^2)+m^8*(8*M^2+5*mchi^2)+m^2*(M^2-mchi^2)^2*(M^4-5*mchi^4)-2*m^6*(8*M^4+9*M^2*mchi^2+5*mchi^4)+2*m^4*(5*M^6+3*M^4*mchi^2+3*M^2*mchi^4+5*mchi^6)-12*m^2*M^4*(m^2-M^2+mchi^2)*Sqrt[m^4-2*m^2*M^2+M^4-2*m^2*mchi^2-2*M^2*mchi^2+mchi^4]*Log[(m^2+M^2-mchi^2+Sqrt[m^4-2*m^2*M^2+M^4-2*m^2*mchi^2-2*M^2*mchi^2+mchi^4])/(2*m*M)]))/(192*M^2*(m^4+(M^2-mchi^2)^2-2*m^2*(M^2+mchi^2))^3*Pi);

gG2[M_,m_,mchi_,Alfas_,gDM_]:=Alfas*gDM^2*mchi^2*(-((m^4+(M^2-mchi^2)^2-2*m^2*(M^2+mchi^2))*(mchi^2*(-2*m^6+(2*M^2-3*mchi^2)*(M^2-mchi^2)^2+m^4*(6*M^2+7*mchi^2)-2*m^2*(3*M^4+mchi^4))+(m^4+(M^2-mchi^2)^2-2*m^2*(M^2+mchi^2))^2*Log[m^2/M^2]))+2*Sqrt[m^4-2*m^2*M^2+M^4-2*m^2*mchi^2-2*M^2*mchi^2+mchi^4]*(m^10-(M^2-mchi^2)^5-5*m^8*(M^2+mchi^2)+10*m^6*(M^4+M^2*mchi^2+mchi^4)-2*m^4*(5*M^6+4*mchi^6)+m^2*(5*M^8-10*M^6*mchi^2+4*M^2*mchi^6+mchi^8))*Log[(m^2+M^2-mchi^2+Sqrt[m^4-2*m^2*M^2+M^4-2*m^2*mchi^2-2*M^2*mchi^2+mchi^4])/(2*m*M)])/(48*mchi^5*(m^4+(M^2-mchi^2)^2-2*m^2*(M^2+mchi^2))^3*Pi);

gG1[M_,m_,mchi_,Alfas_,gDM_]:=Alfas*gDM^2*mchi*(-((m^4+(M^2-mchi^2)^2-2*m^2*(M^2+mchi^2))*(2*mchi^2*(-m^2+M^2+3*mchi^2)+(m^4+(M^2-mchi^2)^2-2*m^2*(M^2+mchi^2))*Log[m^2/M^2]))+2*Sqrt[m^4-2*m^2*M^2+M^4-2*m^2*mchi^2-2*M^2*mchi^2+mchi^4]*(m^6-M^6+3*M^4*mchi^2+M^2*mchi^4-3*mchi^6-3*m^4*(M^2+mchi^2)+m^2*(3*M^4+5*mchi^4))*Log[(m^2+M^2-mchi^2+Sqrt[m^4-2*m^2*M^2+M^4-2*m^2*mchi^2-2*M^2*mchi^2+mchi^4])/(2*m*M)])/(192*mchi^4*(m^4+(M^2-mchi^2)^2-2*m^2*(M^2+mchi^2))^2*Pi);

fQ[M_,m_,mchi_,Alfas_,gDM_]:=gDM^2*mchi/16/(M^2 - (mchi+m)^2)^2;
g1Q[M_,m_,mchi_,Alfas_,gDM_]:=gDM^2*mchi/8/(M^2 - (mchi+m)^2)^2;
g2Q[M_,m_,mchi_,Alfas_,gDM_]:=0;



(* ::Subsubtitle::Closed:: *)
(*Define Nuclear cross - sections*)


(*Define Nuclear cross - sections*)
(*Returns the nucleon matrix elements as a list: {Subscript[f, p],Subscript[f, n]}*)
fnpbymnp[M_,m_,mchi_,Alfas_,gDM_]:=Module[{mli=m,fQ1,fQ2,retp,retn,fG1,fG2},fQ1=Table[fQ[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,1,6}];
fQ2=Table[g1Q[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,1,6}];

fG1=Sum[fG[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,1,6}];
fG2=Sum[gG1[M,mli[[i]],mchi,Alfas,gDM[[i]]] + gG2[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,1,6}];
retp=Sum[fQ1[[i]]*fTqlist[[i]],{i,1,3}]+ 2/27 Sum[fQ1[[i]]*fTG,{i,3,6}]+3/4 *Sum[fQ2[[i]]*Q2listZ[[i]],{i,1,6}]- (8 Pi/ (9* Alfas))fG1*fTG + 3/4 fG2*G2;retn=Sum[fQ1[[i]]*fTqlistn[[i]],{i,1,3}]+ 2/27 Sum[fQ1[[i]]*fTG,{i,3,6}]+3/4 *Sum[fQ2[[i]]*Q2listnZ[[i]],{i,1,6}]- (8 Pi/ (9* Alfas))fG1*fTG + 3/4 fG2*G2;{retp,retn}];

(***************Including RGE**********************)
fnpbymnpwRGE[M_,m_,mchi_,Alfas_,gDM_]:=Module[{fQ1,fQ2,fG1,fG2,l0,l2,retp,retn,mli=m,(*nn=Length[m],*)c0l,c2l},
fG1=Sum[fG[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,1,6}];
fG2=Sum[gG1[M,mli[[i]],mchi,Alfas,gDM[[i]]] + gG2[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,5,6}];

fQ1=Table[fQ[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,1,6}];
fQ2=Table[g1Q[M,mli[[i]],mchi,Alfas,gDM[[i]]],{i,1,6}];
c0l=Table[0,{i,1,7}];
c2l=Table[0,{i,1,7}];
Do[c0l[[i]]=fQ1[[i]];c2l[[i]]=fQ2[[i]],{i,1,6}];
c0l[[7]]=fG1;
c2l[[7]]=fG2;
l0=renc0mu[c0l,2];
l2=renc2mu[c2l,2];
(*Print[l0,"\t",l2];*)

retp=Sum[l0[[i]]*fTqlist[[i]],{i,1,3}] + 2/27 fQ1[[4]]*fTG+ 3/4 *Sum[ l2[[i]]*Q2list[[i]],{i,1,3}]- (8 Pi/ (9* Alfas))l0[[5]]*fTG + 3/4 (l2[[5]])*G2;

retn=Sum[l0[[i]]*fTqlistn[[i]],{i,1,3}] + 2/27 fQ1[[4]]*fTG+ 3/4 *Sum[ l2[[i]]*Q2listn[[i]],{i,1,3}]- (8 Pi/ (9* Alfas))l0[[5]]*fTG + 3/4 (l2[[5]])*G2;{retp,retn}];




(*Define Conversions*)
ConvGEVtoPB=0.389379338* 10^9;
ConvGEVtoCM2=ConvGEVtoPB* 10 ^ -36 ;
(*Define cross-sections*)
sigp[M_,m_,mchi_,Alfas_,gDM_]:=((4/Pi) (mchi*mp/(mchi+mp))^2 * (Abs[fNbymN[M,m,mchi,Alfas,gDM]]*mp)^2)* ConvGEVtoCM2;

sigpRGE[M_,m_,mchi_,Alfas_,gDM_]:=((4/Pi) (mchi*mp/(mchi+mp))^2 * (Abs[fNbymNwRGE[M,m,mchi,Alfas,gDM]]*mp)^2)*ConvGEVtoCM2;

sigSI[M_,m_,mchi_,Alfas_,gDM_]:=Module[{fpn=fnpbymnp[M,m,mchi,Alfas,gDM],mT=mp*Axe},((4/Pi) (mchi*mT/(mchi+mT))^2 * (Abs[Zxe*fpn[[1]]*mp + (Axe-Zxe)*fpn[[2]]*mp ])^2)* ConvGEVtoCM2];

sigSIRGE[M_,m_,mchi_,Alfas_,gDM_]:=Module[{fpn=fnpbymnpwRGE[M,m,mchi,Alfas,gDM],mT=mp*Axe},((4/Pi) (mchi*mT/(mchi+mT))^2 * (Abs[Zxe*fpn[[1]]*mp + (Axe-Zxe)*fpn[[2]]*mp ])^2)* ConvGEVtoCM2];

sigSDuR[M_,m_,mchi_,Alfas_,gDM_]:=3/(64 Pi) (mchi*mp/(mchi+mp))^2 gDM^4 /(M^2 - mchi^2)^2  \[CapitalDelta]u^2;

sigSIuRdirac[M_,m_,mchi_,Alfas_,gDM_]:=1/(64 Pi) (mchi*mp/(mchi+mp))^2 gDM^4 /(M^2 - mchi^2)^2  (1 + Zxe/Axe)^2;

fr3RGESIurdirac[msq_,mchi_]:=Module[{ggdm={1,0,0,0,0,0},tt},If[msq>mchi,tt=(Exp[x1tdt[mchi]]/sigSIuRdirac[msq,mu,mchi,\[Alpha]srun[1.2],1])^(1/4),tt=100];tt];


(* ::Subsubtitle:: *)
(*Import SI limit data *)


(*Import Xenon and Darwin Data*)
SetDirectory[NotebookDirectory[]];
x1tdata=Import["Import/x1t_firstdata_result.csv","CSV"];
x1tdt=Interpolation[Table[{x1tdata[[i,1]],Log[x1tdata[[i,2]]*10^(-45)]},{i,2,Length[x1tdata]}],InterpolationOrder->3,Method-> "Spline"];
x1SIdt=Interpolation[Table[{x1tdata[[i,1]],Log[( ((mchi*mT/(mchi+mT))^2 /  (mchi*mp/(mchi+mp))^2 /.{mT-> mp*Axe,mchi-> x1tdata[[i,1]]}))Axe^2*x1tdata[[i,2]]*10^(-45)]},{i,2,Length[x1tdata]}],InterpolationOrder->3,Method-> "Spline"];
(********)
(*Extrapolate*)
x1dtExtrap=Interpolation[Table[{Log10[x1tdata[[i,1]]],Log[x1tdata[[i,2]]*10^(-45)]},{i,2,Length[x1tdata]}],InterpolationOrder->1,Method-> "Spline"];
x1SIdtExtrap=Interpolation[Table[{Log10[x1tdata[[i,1]]],Log[( ((mchi*mT/(mchi+mT))^2 /  (mchi*mp/(mchi+mp))^2 /.{mT-> mp*Axe,mchi-> x1tdata[[i,1]]}))Axe^2*x1tdata[[i,2]]*10^(-45)]},{i,2,Length[x1tdata]}],InterpolationOrder->1,Method-> "Spline"];



x1dtExtrap2= Interpolation[Table[{1*10^(i),x1dtExtrap[i]},{i,0.1,5,0.05}],InterpolationOrder->1,Method-> "Spline"];

x1SIdtExtrap2= Interpolation[Table[{1*10^(i),x1SIdtExtrap[i]},{i,0.1,5,0.05}],InterpolationOrder->1,Method-> "Spline"];
(********************************************)
(*Import constraints from Darwin in cm^2*)
darwindata=Import["Import/JZ2_Darwin.csv","CSV"];
darwindt=Interpolation[Table[{darwindata[[i,1]],Log[darwindata[[i,2]]]},{i,1,Length[darwindata]}],InterpolationOrder->3,Method-> "Spline"];
darwinSIdt=Interpolation[Table[{darwindata[[i,1]],Log[Axe^2*darwindata[[i,2]]]},{i,1,Length[darwindata]}],InterpolationOrder->3,Method-> "Spline"];

darwindtExtrap=Interpolation[Table[{Log10[darwindata[[i,1]]],Log[darwindata[[i,2]]]},{i,1,Length[darwindata]}],InterpolationOrder->1,Method-> "Spline"];
darwinSIdtExtrap=Interpolation[Table[{Log10[darwindata[[i,1]]],Log[( ((mchi*mT/(mchi+mT))^2 /  (mchi*mp/(mchi+mp))^2 /.{mT-> mp*Axe,mchi-> darwindata[[i,1]]}))Axe^2*darwindata[[i,2]]]},{i,1,Length[darwindata]}],InterpolationOrder->1,Method-> "Spline"];


darwindtExtrap2= Interpolation[Table[{1*10^(i),darwindtExtrap[i]},{i,0.1,5,0.05}],InterpolationOrder->1,Method-> "Spline"];

darwinSIdtExtrap2= Interpolation[Table[{1*10^(i),darwinSIdtExtrap[i]},{i,0.1,5,0.05}],InterpolationOrder->1,Method-> "Spline"];


(* ::Subsubtitle:: *)
(*Import SD limit data *)


SetDirectory[NotebookDirectory[]];
sdtdata=Import["Import/pico60.dat"];
sdtdt=Interpolation[Table[{sdtdata[[i,1]],Log[sdtdata[[i,2]]]},{i,2,Length[sdtdata]}],InterpolationOrder->3,Method-> "Spline"];

DARWINsdtdata=Import["Import/DARWIN_SD.dat"];
DARWINsdtdt=Interpolation[Table[{sdtdata[[i,1]],Log[sdtdata[[i,2]]*10^(36)]},{i,2,Length[sdtdata]}],InterpolationOrder->3,Method-> "Spline"];

DarwinSDdtExtrap=Interpolation[Table[{Log10[DARWINsdtdata[[i,1]]],Log[DARWINsdtdata[[i,2]]*10^(36)]},{i,1,Length[DARWINsdtdata]}],InterpolationOrder->1,Method-> "Spline"];
DarwinSDdtExtrap2= Interpolation[Table[{1*10^(i),DarwinSDdtExtrap[i]},{i,0.1,5,0.05}],InterpolationOrder->1,Method-> "Spline"];



(* ::Subsubtitle:: *)
(*Define SD nucleon cross-section with RGE*)


RASinglet[nf_,\[Mu]l_,\[Mu]h_]:=Module[{\[Beta]0=11 -2/3nf},Exp[2 nf /(Pi \[Beta]0)(\[Alpha]srun[\[Mu]h]-\[Alpha]srun[\[Mu]l])]];
(*Input gDM as a list of length 3 (u,d,s). 
*)
sdwilly[M_,mchi_,g_]:=g^2/(M^2 - mchi^2);
sigSDp[M_,mchi_,gDM_]:=Module[{cA0,cA3,cA8,cu,cd,cs,sigp},

 cu=sdwilly[M[[1]],mchi,gDM[[1]]];cd=sdwilly[M[[2]],mchi,gDM[[2]]];cs=sdwilly[M[[3]],mchi,gDM[[3]]];

cA0=cu+cd+cs;cA3=cu-cd;cA8=1/(Sqrt[3])(cu+cd -2 cs);

sigp=3/(16 Pi) (mchi*mp/(mchi+mp))^2 (cu \[CapitalDelta]up + cd  \[CapitalDelta]dp + cs \[CapitalDelta]sp)^2;

(*Print["cu= ",cu,"\t cd=",cd,"\t cs=",cs," sigp=" ,sigp];
Print["cA0=",cA0," cA3=",cA3," cA8=",cA8];*)

cA0=RASinglet[3,mul,mc]RASinglet[4,mc,mb] RASinglet[5,mb,mt]cA0;(*Print["cA0R=",cA0] ;*)cu=1/6 (2 cA0+3 cA3+Sqrt[3] cA8) ;cd=1/6 (2 cA0-3 cA3+Sqrt[3] cA8);cs=1/3 (cA0-Sqrt[3] cA8);sigp=3/(16 Pi) (mchi*mp/(mchi+mp))^2 (cu \[CapitalDelta]up + cd  \[CapitalDelta]dp + cs \[CapitalDelta]sp)^2;(*Print["cu= ",cu,"\t cd=",cd,"\t cs=",cs," sigp=" ,sigp];*)sigp*ConvGEVtoPB ];


(********SD cross-sections*********)
sigSDuR[M_,m_,mchi_,Alfas_,gDM_,\[CapitalDelta]u_]:=3/(16 Pi) (mchi*mp/(mchi+mp))^2 gDM^4 /(M^2 - mchi^2)^2  \[CapitalDelta]u^2 * ConvGEVtoPB;

sigSDdR[M_,m_,mchi_,Alfas_,gDM_,\[CapitalDelta]d_,\[CapitalDelta]s_]:=3/(16 Pi) (mchi*mp/(mchi+mp))^2 gDM^4 /(M^2 - mchi^2)^2  (\[CapitalDelta]d + \[CapitalDelta]s)^2*ConvGEVtoPB;

sigSDqL[M_,m_,mchi_,Alfas_,gDM_,\[CapitalDelta]u_,\[CapitalDelta]d_,\[CapitalDelta]s_]:=3/(16 Pi) (mchi*mp/(mchi+mp))^2 gDM^4 /(M^2 - mchi^2)^2  (\[CapitalDelta]u+\[CapitalDelta]d + \[CapitalDelta]s)^2* ConvGEVtoPB;

(******WITHOUT RGE*****)

fr3SDuR[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[sdtdt[mchi]]/sigSDuR[msq,mu,mchi,\[Alpha]srun[mul],1,\[CapitalDelta]up])^(1/4),tt=100];tt];

fr3SDdR[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[sdtdt[mchi]]/sigSDdR[msq,mu,mchi,\[Alpha]srun[mul],1,\[CapitalDelta]dp,\[CapitalDelta]sp])^(1/4),tt=100];tt];

fr3SDqL[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[sdtdt[mchi]]/sigSDqL[msq,mu,mchi,\[Alpha]srun[mul],1,\[CapitalDelta]up,\[CapitalDelta]dp,\[CapitalDelta]sp])^(1/4),tt=100];tt];
(*****WITH RGE*******)
fr3RGESDuR[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[sdtdt[mchi]]/sigSDp[msq{1,1,1},mchi,{1,0,0}])^(1/4),tt=100];tt];

fr3RGESDdR[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[sdtdt[mchi]]/sigSDp[msq{1,1,1},mchi,{0,1,1}])^(1/4),tt=100];tt];

fr3RGESDqL[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[sdtdt[mchi]]/sigSDp[msq{1,1,1},mchi,{1,1,1}])^(1/4),tt=100];tt];


(*****WITH RGE*******)
Darwinfr3RGESDuR[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[DarwinSDdtExtrap2[mchi]]/sigSDp[msq{1,1,1},mchi,{1,0,0}])^(1/4),tt=100];tt];

Darwinfr3RGESDdR[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[DarwinSDdtExtrap2[mchi]]/sigSDp[msq{1,1,1},mchi,{0,1,1}])^(1/4),tt=100];tt];

Darwinfr3RGESDqL[msq_?NumericQ,mchi_?NumericQ]:=Module[{tt=0},If[msq>mchi,tt=(Exp[DarwinSDdtExtrap2[mchi]]/sigSDp[msq{1,1,1},mchi,{1,1,1}])^(1/4),tt=100];tt];
