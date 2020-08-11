function f = SRMQ1C2cfa(IMAGE,Tc,part)
% -------------------------------------------------------------------------
% Copyright (c) 2014 DDE Lab, Binghamton University, NY.
% All Rights Reserved.
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The program is supplied "as is," without any
% accompanying services from DDE Lab. DDE Lab does not warrant the
% operation of the program will be uninterrupted or error-free. The
% end-user understands that the program was developed for research purposes
% and is advised not to rely exclusively on the program for any reason. In
% no event shall Binghamton University or DDE Lab be liable to any party
% for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software. DDE Lab
% disclaims any warranties, and has no obligations to provide maintenance,
% support, updates, enhancements or modifications.
% -------------------------------------------------------------------------
% Contact: mgoljan@binghamton.edu | fridrich@binghamton.edu | May 2014
%          http://dde.binghamton.edu/download/feature_extractors
% -------------------------------------------------------------------------
% SRMQ1C2cfa is an extension of SRMQ1 for color images with CFA artifacts 
%   All residuals are split to green and red/blue pixels. 
%   Features dimension thus doubles.
% -------------------------------------------------------------------------
% Extracts all 106 submodels presented in [1] as part of a rich model for
% steganalysis of digital images. All features are calculated in the
% spatial domain, summed over three color channels. 3D inter-channel 
% cooccurrences are computed for each residual and feature vectors are 
% stored in a structured variable 'f'. 
% -------------------------------------------------------------------------
% Input:  IMAGE ... path to the image (can be JPEG)
%         Tc ...... Truncation for quantized residuals that enter CRMQ1
%                   Truncation in SRMQ1 is fixed to T=2
%         part .... 'color' or anything other than 'all' -output only CRMQ1
%                   otherwise output all SCRMQ1 features
% Output: f ....... extracted SRM features in a structured format
% -------------------------------------------------------------------------
% [1] Rich Models for Steganalysis of Digital Images, J. Fridrich and J.
% Kodovsky, IEEE Transactions on Information Forensics and Security, 2011.
% [2] CFA-aware Features for Steganalysis of Color Images, M. Goljan and 
% J. Fridrich, Proceedings SPIE, Electronic Imaging, Media Forensics
% and Security, 2015
% -------------------------------------------------------------------------
if ischar(IMAGE)
    X = double(imread(IMAGE));
else
    X = double(IMAGE);
end
f = [];
if nargin<3,
    part = 'all';
elseif isempty(part),
    part = 'all';
end
if nargin<2,
    % SRMC for color image, T=3 fixed
    ,Tc=3,      % truncation for between color cooccurrences
end
f = post_processing(all1st3rdc(X,1,Tc,'1st'),'f1',1,f,Tc);   % 1st order, q=1
% f = post_processing(all1stc(X,2,Tc),'f1',2,f,Tc); % 1st order, q=2
for q=[1], f = post_processing(all2ndc(X,q*2,Tc),'f2',q,f,Tc); end    % 2nd order
for q=[1], f = post_processing(all1st3rdc(X,q*3,Tc,'3rd'),'f3',q,f,Tc); end 	% 3rd order
for q=[1], f = post_processing(all3x3c(X,q*4,Tc),'f3x3',q,f,Tc); end  % 3x3
for q=[1], f = post_processing(all5x5c(X,q*12,Tc),'f5x5',q,f,Tc); end % 5x5

if strcmp(part,'all')       % do all submodels
% SRM or merged SRM for color image, T=2 fixed
f = post_processing(all1st3rd(X,1,'1st'),'f1',1,f);   % 1st order, q=1
% f = post_processing(all1st(X,2),'f1',2,f); % 1st order, q=2
for q=[1], f = post_processing(all2nd(X,q*2),'f2',q,f); end    % 2nd order
for q=[1], f = post_processing(all1st3rd(X,q*3,'3rd'),'f3',q,f); end 	% 3rd order
for q=[1], f = post_processing(all3x3(X,q*4),'f3x3',q,f); end  % 3x3
for q=[1], f = post_processing(all5x5(X,q*12),'f5x5',q,f); end % 5x5
end                         % do all submodels

%**** FUNCTION ****%
function RESULT = post_processing(DATA,f,q,RESULT,Tc)

Ss = fieldnames(DATA);
for Sid = 1:length(Ss)
    VARNAME = [f '_' Ss{Sid} '_q' strrep(num2str(q),'.','')];
    eval(['RESULT.' VARNAME ' = reshape(single(DATA.' Ss{Sid} '),1,[]);' ])
end

% symmetrize
L = fieldnames(RESULT);
for i=1:length(L)
    name = L{i}; % feature name
    if name(1)=='s', continue; end
    [T,N,Q] = parse_feaname(name);
    if strcmp(T,''), continue; end
    % symmetrization
    if strcmp(N(1:3),'min') || strcmp(N(1:3),'max')
        % minmax symmetrization
        OUT = ['s' T(2:end) '_minmax' N(4:end) '_' Q];
        if isfield(RESULT,OUT), continue; end
        eval(['Fmin = RESULT.' strrep(name,'max','min') ';']);
        eval(['Fmax = RESULT.' strrep(name,'min','max') ';']);
        if strcmp(N(end),'c'),       	% color minmax symmetrization, T=3
%             F = symfea([Fmin Fmax]',Tc,3,'mnmx')';   % including rgb<->bgr symm
            F = Symm2split([Fmin Fmax],Tc,3,'minmaxdir');
        else
%             F = symfea([Fmin Fmax]',2,4,'mnmx')';   %#ok<*NASGU>
            F = Symm2split([Fmin Fmax],2,4,'minmaxdir');
        end
        eval(['RESULT.' OUT ' = single(F);' ]);
    elseif strcmp(N(1:4),'spam')
        % spam symmetrization
        OUT = ['s' T(2:end) '_' N '_' Q];
        if isfield(RESULT,OUT), continue; end
        eval(['Fold = RESULT.' name ';']);
        if strcmp(N(end),'c'),     	   % color sign symmetrization, T=3
%             F = symm1sign(Fold',3,3)'; % w/o rgb<->bgr symm
%             F = symm1(Fold',Tc,3)';    % including rgb<->bgr symm
            F = Symm2split(Fold,Tc,3,'dirsign');    % including rgb<->bgr symm
        else
%             F = symm1(Fold',2,4)';
            F = Symm2split(Fold,2,4,'dirsign');
        end
        eval(['RESULT.' OUT ' = single(F);' ]);
    end
end


% delete RESULT.f*
L = fieldnames(RESULT);
for i=1:length(L)
    name = L{i}; % feature name
    if name(1)=='f'
        RESULT = rmfield(RESULT,name);
    end
end
% merge spam features
L = fieldnames(RESULT);
for i=1:length(L)
    name = L{i}; % feature name
    [T,N,Q] = parse_feaname(name);
    if ~strcmp(N(1:4),'spam'), continue; end
    if strcmp(T,''), continue; end
    if strcmp(N(end),'v')||(strcmp(N,'spam11')&&strcmp(T,'s5x5'))
    elseif strcmp(N(end),'h')
        % h+v union
        OUT = [T '_' N 'v_' Q ];
        if isfield(RESULT,OUT), continue; end
        name2 = strrep(name,'h_','v_');
        eval(['Fh = RESULT.' name ';']);
        eval(['Fv = RESULT.' name2 ';']);
        eval(['RESULT.' OUT ' = [Fh Fv];']);
        RESULT = rmfield(RESULT,name);
        RESULT = rmfield(RESULT,name2);
    elseif strcmp(N,'spam11')||strcmp(N,'spam11c')
        % KBKV creation
        OUT = ['s35_' N '_' Q];
        if isfield(RESULT,OUT), continue; end
        name1 = strrep(name,'5x5','3x3');
        name2 = strrep(name,'3x3','5x5');
        if ~isfield(RESULT,name1), continue; end
        if ~isfield(RESULT,name2), continue; end
        eval(['F_KB = RESULT.' name1 ';']);
        eval(['F_KV = RESULT.' name2 ';']);
        eval(['RESULT.' OUT ' = [F_KB F_KV];']);
        RESULT = rmfield(RESULT,name1);
        RESULT = rmfield(RESULT,name2);
    end
end

%**** FUNCTION ****%
function [T,N,Q] = parse_feaname(name)
[T,N,Q] = deal('');
P = strfind(name,'_'); if length(P)~=2, return; end
T = name(1:P(1)-1);
N = name(P(1)+1:P(2)-1);
Q = name(P(2)+1:end);

%**** FUNCTION ****%
function g = all1st3rd(X,q,resorder)
%
% SRM coccurrence features of 1st or 3rd order
% resorder  = '1st' or '3rd'
% X must be a matrix of doubles or singles (the image) and q is the 
% quantization step (any positive number).
%
% Recommended values of q are c, 1.5c, 2c, where c is the central
% coefficient in the differential (at X(I,J)).
%
% This function outputs co-occurrences of ALL 1st-order residuals
% listed in Figure 1 in our journal HUGO paper (version from June 14), 
% including the naming convention.
%
% List of outputted features:
%
% 1a) spam14h
% 1b) spam14v (orthogonal-spam)
% 1c) minmax22v
% 1d) minmax24
% 1e) minmax34v
% 1f) minmax41
% 1g) minmax34
% 1h) minmax48h
% 1i) minmax54
%
% Naming convention:
%
% name = {type}{f}{sigma}{scan}
% type \in {spam, minmax}
% f \in {1,2,3,4,5} number of filters that are "minmaxed"
% sigma \in {1,2,3,4,8} symmetry index
% scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
% h and v scans).
%
% All odd residuals are implemented the same way simply by
% narrowing the range for I and J and replacing the residuals --
% -- they should "stick out" (trcet) in the same direction as 
% the 1st order ones. For example, for the 3rd order:
%
% RU = -X(I-2,J+2)+3*X(I-1,J+1)-3*X(I,J)+X(I+1,J-1); ... etc.
%
% Note1: The term X(I,J) should always have the "-" sign.
% Note2: This script does not include s, so, cout, cin versions (weak).
% This function calls Cooc1split.m, Cooc1.m, Quant.m, and more

[M N three] = size(X);  [T,order] = deal(2,4);
    switch resorder
        case '1st'
         [I,J] = deal(2:M-1,2:N-1); 
        case '3rd'
         [I,J] = deal(3:M-2,3:N-2); 
    end

% Variable names are self-explanatory (R = right, U = up, L = left, D = down)
for j=1:three
    switch resorder
        case '1st'
         [R,L,U,D] = deal(X(I,J+1,j)-X(I,J,j),X(I,J-1,j)-X(I,J,j),X(I-1,J,j)-X(I,J,j),X(I+1,J,j)-X(I,J,j)); 
         [RU{j},LU{j},RD{j},LD{j}] = deal(X(I-1,J+1,j)-X(I,J,j),X(I-1,J-1,j)-X(I,J,j),X(I+1,J+1,j)-X(I,J,j),X(I+1,J-1,j)-X(I,J,j));
        case '3rd'
         [R,L,U,D] = deal(-X(I,J+2,j)+3*X(I,J+1,j)-3*X(I,J,j)+X(I,J-1,j),-X(I,J-2,j)+3*X(I,J-1,j)-3*X(I,J,j)+...
            X(I,J+1,j),-X(I-2,J,j)+3*X(I-1,J,j)-3*X(I,J,j)+X(I+1,J,j),-X(I+2,J,j)+3*X(I+1,J,j)-3*X(I,J,j)+X(I-1,J,j));
         [RU{j},LU{j},RD{j},LD{j}] = deal(-X(I-2,J+2,j)+3*X(I-1,J+1,j)-3*X(I,J,j)+X(I+1,J-1,j),-X(I-2,J-2,j)+3*X(I-1,J-1,j)-3*X(I,J,j)+...
            X(I+1,J+1,j),-X(I+2,J+2,j)+3*X(I+1,J+1,j)-3*X(I,J,j)+X(I-1,J-1,j),-X(I+2,J-2,j)+3*X(I+1,J-1,j)-3*X(I,J,j)+X(I-1,J+1,j));
    end
    [Rq{j},Lq{j},Uq{j},Dq{j}] = deal(Quant(R,q,T),Quant(L,q,T),Quant(U,q,T),Quant(D,q,T));
    [RU{j},LU{j},RD{j},LD{j}] = deal(Quant(RU{j},q,T),Quant(LU{j},q,T),Quant(RD{j},q,T),Quant(LD{j},q,T));
end,       clear R L U D X

% spam14h/v -- to be symmetrized as spam, directional, hv-nonsymmetrical
 g.spam14h = Cooc1split(Rq,order,'hor',T) + Cooc1split(Uq,order,'ver',T);
 g.spam14v = Cooc1split(Rq,order,'ver',T) + Cooc1split(Uq,order,'hor',T);

% minmax22h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical.
for j=1:three
 [RL_min{j},UD_min{j},RL_max{j},UD_max{j}] = deal(min(Rq{j},Lq{j}),min(Uq{j},Dq{j}),max(Rq{j},Lq{j}),max(Uq{j},Dq{j}));
end
 g.min22h = Cooc1split(RL_min,order,'hor',T) + Cooc1split(UD_min,order,'ver',T);
 g.max22h = Cooc1split(RL_max,order,'hor',T) + Cooc1split(UD_max,order,'ver',T);

% minmax22v -- to be symmetrized as mnmx, directional, hv-nonsymmetrical. Good with higher-order residuals! Note: 22h is bad (too much neighborhood overlap).
 g.min22v = Cooc1split(RL_min,order,'ver',T) + Cooc1split(UD_min,order,'hor',T);
 g.max22v = Cooc1split(RL_max,order,'ver',T) + Cooc1split(UD_max,order,'hor',T);

% minmax24 -- to be symmetrized as mnmx, directional, hv-symmetrical. Darn good, too.
for j=1:three
 [RU_min{j},RD_min{j},LU_min{j},LD_min{j}] = deal(min(Rq{j},Uq{j}),min(Rq{j},Dq{j}),min(Lq{j},Uq{j}),min(Lq{j},Dq{j}));
 [RU_max{j},RD_max{j},LU_max{j},LD_max{j}] = deal(max(Rq{j},Uq{j}),max(Rq{j},Dq{j}),max(Lq{j},Uq{j}),max(Lq{j},Dq{j}));
end
 g.min24 = Cooc1split(vercat(vercat(RU_min,RD_min),vercat(LU_min,LD_min)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RU_min,RD_min),horcat(LU_min,LD_min)),order,'ver',T);
 g.max24 = Cooc1split(vercat(vercat(RU_max,RD_max),vercat(LU_max,LD_max)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RU_max,RD_max),horcat(LU_max,LD_max)),order,'ver',T);

% minmax34h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical
for j=1:three
 [Uq_min{j},Rq_min{j},Dq_min{j},Lq_min{j}] = deal(min(min(Lq{j},Uq{j}),Rq{j}),min(min(Uq{j},Rq{j}),Dq{j}),min(min(Rq{j},Dq{j}),Lq{j}),min(min(Dq{j},Lq{j}),Uq{j}));
 [Uq_max{j},Rq_max{j},Dq_max{j},Lq_max{j}] = deal(max(max(Lq{j},Uq{j}),Rq{j}),max(max(Uq{j},Rq{j}),Dq{j}),max(max(Rq{j},Dq{j}),Lq{j}),max(max(Dq{j},Lq{j}),Uq{j}));
end,    clear Rq Uq Dq Lq 
 g.min34h = Cooc1split(vercat(Uq_min,Dq_min),order,'hor',T) + Cooc1split(horcat(Lq_min,Rq_min),order,'ver',T);
 g.max34h = Cooc1split(vercat(Uq_max,Dq_max),order,'hor',T) + Cooc1split(horcat(Rq_max,Lq_max),order,'ver',T);

% minmax34v -- v works well, h does not, to be symmetrized as mnmx, directional, hv-nonsymmetrical
 g.min34v = Cooc1split(horcat(Uq_min,Dq_min),order,'ver',T) + Cooc1split(vercat(Rq_min,Lq_min),order,'hor',T);
 g.max34v = Cooc1split(horcat(Uq_max,Dq_max),order,'ver',T) + Cooc1split(vercat(Rq_max,Lq_max),order,'hor',T);
    clear Uq_min Dq_min Rq_min Lq_min 
    clear Uq_max Dq_max Rq_max Lq_max
 
% minmax41 -- to be symmetrized as mnmx, non-directional, hv-symmetrical
for j=1:three
 [R_min{j},R_max{j}] = deal(min(RL_min{j},UD_min{j}),max(RL_max{j},UD_max{j}));  
end,   	clear RL_min UD_min RL_max UD_max
 g.min41 = Cooc1split(R_min,order,'hor',T) + Cooc1split(R_min,order,'ver',T);
 g.max41 = Cooc1split(R_max,order,'hor',T) + Cooc1split(R_max,order,'ver',T);
    clear R_min R_max

% minmax34 -- good, to be symmetrized as mnmx, directional, hv-symmetrical
for j=1:three
 [RU_min{j},RD_min{j},LU_min{j},LD_min{j}] = deal(min(RU_min{j},RU{j}),min(RD_min{j},RD{j}),min(LU_min{j},LU{j}),min(LD_min{j},LD{j}));
 [RU_max{j},RD_max{j},LU_max{j},LD_max{j}] = deal(max(RU_max{j},RU{j}),max(RD_max{j},RD{j}),max(LU_max{j},LU{j}),max(LD_max{j},LD{j}));
end
 g.min34 = Cooc1split(vercat(vercat(RU_min,RD_min),vercat(LU_min,LD_min)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RU_min,RD_min),horcat(LU_min,LD_min)),order,'ver',T);
 g.max34 = Cooc1split(vercat(vercat(RU_max,RD_max),vercat(LU_max,LD_max)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RU_max,RD_max),horcat(LU_max,LD_max)),order,'ver',T);

 % minmax48h -- h better than v, to be symmetrized as mnmx, directional, hv-nonsymmetrical. 48v is almost as good as 48h; for 3rd-order but weaker for 1st-order. Here, I am outputting both but Figure 1 in our paper lists only 48h.
for j=1:three
 [RU_min2{j},RD_min2{j},LD_min2{j},LU_min2{j}] = deal(min(RU_min{j},LU{j}),min(RD_min{j},RU{j}),min(LD_min{j},RD{j}),min(LU_min{j},LD{j}));
 [RU_min3{j},RD_min3{j},LD_min3{j},LU_min3{j}] = deal(min(RU_min{j},RD{j}),min(RD_min{j},LD{j}),min(LD_min{j},LU{j}),min(LU_min{j},RU{j}));
 [RU_max2{j},RD_max2{j},LD_max2{j},LU_max2{j}] = deal(max(RU_max{j},LU{j}),max(RD_max{j},RU{j}),max(LD_max{j},RD{j}),max(LU_max{j},LD{j}));
 [RU_max3{j},RD_max3{j},LD_max3{j},LU_max3{j}] = deal(max(RU_max{j},RD{j}),max(RD_max{j},LD{j}),max(LD_max{j},LU{j}),max(LU_max{j},RU{j}));
end,    clear RU_min RD_min LU_min LD_min RU_max RD_max LU_max LD_max
 g.min48h = Cooc1split(vercat(vercat(RU_min2,LD_min2),vercat(RD_min3,LU_min3)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RD_min2,LU_min2),horcat(RU_min3,LD_min3)),order,'ver',T);
 g.min48v = Cooc1split(vercat(vercat(RD_min2,LU_min2),vercat(RU_min3,LD_min3)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RU_min2,LD_min2),horcat(RD_min3,LU_min3)),order,'ver',T);
 g.max48h = Cooc1split(vercat(vercat(RU_max2,LD_max2),vercat(RD_max3,LU_max3)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RD_max2,LU_max2),horcat(RU_max3,LD_max3)),order,'ver',T);
 g.max48v = Cooc1split(vercat(vercat(RD_max2,LU_max2),vercat(RU_max3,LD_max3)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RU_max2,LD_max2),horcat(RD_max3,LU_max3)),order,'ver',T);

 % minmax54 -- to be symmetrized as mnmx, directional, hv-symmetrical
for j=1:three
 [RU_min4{j},RD_min4{j},LD_min4{j},LU_min4{j}] = deal(min(RU_min2{j},RD{j}),min(RD_min2{j},LD{j}),min(LD_min2{j},LU{j}),min(LU_min2{j},RU{j}));
 [RU_min5{j},RD_min5{j},LD_min5{j},LU_min5{j}] = deal(min(RU_min3{j},LU{j}),min(RD_min3{j},RU{j}),min(LD_min3{j},RD{j}),min(LU_min3{j},LD{j}));
 [RU_max4{j},RD_max4{j},LD_max4{j},LU_max4{j}] = deal(max(RU_max2{j},RD{j}),max(RD_max2{j},LD{j}),max(LD_max2{j},LU{j}),max(LU_max2{j},RU{j}));
 [RU_max5{j},RD_max5{j},LD_max5{j},LU_max5{j}] = deal(max(RU_max3{j},LU{j}),max(RD_max3{j},RU{j}),max(LD_max3{j},RD{j}),max(LU_max3{j},LD{j}));
end,    clear RU_min2 RD_min2 LD_min2 LU_min2 RU_min3 RD_min3 LD_min3 LU_min3 
        clear RU_max2 RD_max2 LD_max2 LU_max2 RU_max3 RD_max3 LD_max3 LU_max3 
 g.min54 = Cooc1split(vercat(vercat(RU_min4,LD_min4),vercat(RD_min5,LU_min5)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RD_min4,LU_min4),horcat(RU_min5,LD_min5)),order,'ver',T);  
 g.max54 = Cooc1split(vercat(vercat(RU_max4,LD_max4),vercat(RD_max5,LU_max5)),order,'hor',T) + ...
     Cooc1split(horcat(horcat(RD_max4,LU_max4),horcat(RU_max5,LD_max5)),order,'ver',T);
 
%**** FUNCTION ****%
function g = all1st3rdc(X,q,T,resorder)
%
% SRMC coccurrence features of 1st or 3rd order
% resorder = '1st' or '3rd'
%
% X must be a matrix of doubles or singles (the color image) and q is the 
% quantization step (any positive number).
%
% Recommended values of q are c, 1.5c, 2c, where c is the central
% coefficient in the differential (at X(I,J)).
%
% This function outputs co-occurrences of ALL 1st-order residuals
% listed in Figure 1 in our journal HUGO paper (version from June 14), 
% including the naming convention.
%
% List of outputted features:
%
% 1a) spam14h
% 1b) spam14v (orthogonal-spam)
% 1c) minmax22v
% 1d) minmax24
% 1e) minmax34v
% 1f) minmax41
% 1g) minmax34
% 1h) minmax48h
% 1i) minmax54
%
% Naming convention:
%
% name = {type}{f}{sigma}{scan}
% type \in {spam, minmax}
% f \in {1,2,3,4,5} number of filters that are "minmaxed"
% sigma \in {1,2,3,4,8} symmetry index
% scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
% h and v scans).
%
% All odd residuals are implemented the same way simply by
% narrowing the range for I and J and replacing the residuals --
% -- they should "stick out" (trcet) in the same direction as 
% the 1st order ones. For example, for the 3rd order:
%
% RU = -X(I-2,J+2)+3*X(I-1,J+1)-3*X(I,J)+X(I+1,J-1); ... etc.
%
% Note1: The term X(I,J) should always have the "-" sign.
% Note2: This script does not include s, so, cout, cin versions (weak).
% This function calls Cooc1split.m and Quant.m

[M N three] = size(X);  order = 4;
    switch resorder
        case '1st'
         [I,J] = deal(2:M-1,2:N-1);
        case '3rd'
         [I,J] = deal(3:M-2,3:N-2);
    end

% Variable names are self-explanatory (R = right, U = up, L = left, D = down)
for j=1:three
    switch resorder
        case '1st'
         [R,L,U,D] = deal(X(I,J+1,j)-X(I,J,j),X(I,J-1,j)-X(I,J,j),X(I-1,J,j)-X(I,J,j),X(I+1,J,j)-X(I,J,j)); 
         [RU{j},LU{j},RD{j},LD{j}] = deal(X(I-1,J+1,j)-X(I,J,j),X(I-1,J-1,j)-X(I,J,j),X(I+1,J+1,j)-X(I,J,j),X(I+1,J-1,j)-X(I,J,j));
        case '3rd'
         [R,L,U,D] = deal(-X(I,J+2,j)+3*X(I,J+1,j)-3*X(I,J,j)+X(I,J-1,j),-X(I,J-2,j)+3*X(I,J-1,j)-3*X(I,J,j)+...
            X(I,J+1,j),-X(I-2,J,j)+3*X(I-1,J,j)-3*X(I,J,j)+X(I+1,J,j),-X(I+2,J,j)+3*X(I+1,J,j)-3*X(I,J,j)+X(I-1,J,j));
         [RU{j},LU{j},RD{j},LD{j}] = deal(-X(I-2,J+2,j)+3*X(I-1,J+1,j)-3*X(I,J,j)+X(I+1,J-1,j),-X(I-2,J-2,j)+3*X(I-1,J-1,j)-3*X(I,J,j)+...
            X(I+1,J+1,j),-X(I+2,J+2,j)+3*X(I+1,J+1,j)-3*X(I,J,j)+X(I-1,J-1,j),-X(I+2,J-2,j)+3*X(I+1,J-1,j)-3*X(I,J,j)+X(I-1,J+1,j));
    end
    [Rq{j},Lq{j},Uq{j},Dq{j}] = deal(Quant(R,q,T),Quant(L,q,T),Quant(U,q,T),Quant(D,q,T));
    [RU{j},LU{j},RD{j},LD{j}] = deal(Quant(RU{j},q,T),Quant(LU{j},q,T),Quant(RD{j},q,T),Quant(LD{j},q,T));
end,       clear R L U D X

% spam14h/v -- to be symmetrized as spam, directional, hv-nonsymmetrical
 g.spam14c = Cooc1split(Rq,3,'col',T) + Cooc1split(Uq,3,'col',T);

% minmax22h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical.
for j=1:three
 [RL_min{j},UD_min{j},RL_max{j},UD_max{j}] = deal(min(Rq{j},Lq{j}),min(Uq{j},Dq{j}),max(Rq{j},Lq{j}),max(Uq{j},Dq{j}));
end
 g.min22c = Cooc1split(RL_min,3,'col',T) + Cooc1split(UD_min,3,'col',T);
 g.max22c = Cooc1split(RL_max,3,'col',T) + Cooc1split(UD_max,3,'col',T);

% minmax24 -- to be symmetrized as mnmx, directional, hv-symmetrical. Darn good, too.
for j=1:three
 [RU_min{j},RD_min{j},LU_min{j},LD_min{j}] = deal(min(Rq{j},Uq{j}),min(Rq{j},Dq{j}),min(Lq{j},Uq{j}),min(Lq{j},Dq{j}));
 [RU_max{j},RD_max{j},LU_max{j},LD_max{j}] = deal(max(Rq{j},Uq{j}),max(Rq{j},Dq{j}),max(Lq{j},Uq{j}),max(Lq{j},Dq{j}));
end
 g.min24c = Cooc1split(vercat(vercat(RU_min,RD_min),vercat(LU_min,LD_min)),3,'col',T);
 g.max24c = Cooc1split(vercat(vercat(RU_max,RD_max),vercat(LU_max,LD_max)),3,'col',T);

% minmax34h -- to be symmetrized as mnmx, directional, hv-nonsymmetrical
for j=1:three
 [Uq_min{j},Rq_min{j},Dq_min{j},Lq_min{j}] = deal(min(min(Lq{j},Uq{j}),Rq{j}),min(min(Uq{j},Rq{j}),Dq{j}),min(min(Rq{j},Dq{j}),Lq{j}),min(min(Dq{j},Lq{j}),Uq{j}));
 [Uq_max{j},Rq_max{j},Dq_max{j},Lq_max{j}] = deal(max(max(Lq{j},Uq{j}),Rq{j}),max(max(Uq{j},Rq{j}),Dq{j}),max(max(Rq{j},Dq{j}),Lq{j}),max(max(Dq{j},Lq{j}),Uq{j}));
end,    clear Rq Uq Dq Lq 

% minmax34v -- v works well, h does not, to be symmetrized as mnmx, directional, hv-nonsymmetrical
 g.min34hvc = Cooc1split(horcat(Uq_min,Dq_min),3,'col',T) + Cooc1split(vercat(Rq_min,Lq_min),3,'col',T);
    clear Uq_min Dq_min Rq_min Lq_min 
 g.max34hvc = Cooc1split(horcat(Uq_max,Dq_max),3,'col',T) + Cooc1split(vercat(Rq_max,Lq_max),3,'col',T);
    clear Uq_max Dq_max Rq_max Lq_max
 
% minmax41 -- to be symmetrized as mnmx, non-directional, hv-symmetrical
for j=1:three
 [R_min{j},R_max{j}] = deal(min(RL_min{j},UD_min{j}),max(RL_max{j},UD_max{j}));  
end,   	clear RL_min UD_min RL_max UD_max
 g.min41c = Cooc1split(R_min,3,'col',T);      clear R_min
 g.max41c = Cooc1split(R_max,3,'col',T);      clear R_max

% minmax34 -- good, to be symmetrized as mnmx, directional, hv-symmetrical
for j=1:three
 [RU_min{j},RD_min{j},LU_min{j},LD_min{j}] = deal(min(RU_min{j},RU{j}),min(RD_min{j},RD{j}),min(LU_min{j},LU{j}),min(LD_min{j},LD{j}));
 [RU_max{j},RD_max{j},LU_max{j},LD_max{j}] = deal(max(RU_max{j},RU{j}),max(RD_max{j},RD{j}),max(LU_max{j},LU{j}),max(LD_max{j},LD{j}));
end
 g.min34c = Cooc1split(vercat(vercat(RU_min,RD_min),vercat(LU_min,LD_min)),3,'col',T);
 g.max34c = Cooc1split(vercat(vercat(RU_max,RD_max),vercat(LU_max,LD_max)),3,'col',T);

 % minmax48h -- h better than v, to be symmetrized as mnmx, directional, hv-nonsymmetrical. 48v is almost as good as 48h; for 3rd-order but weaker for 1st-order. Here, I am outputting both but Figure 1 in our paper lists only 48h.
for j=1:three
 [RU_min2{j},RD_min2{j},LD_min2{j},LU_min2{j}] = deal(min(RU_min{j},LU{j}),min(RD_min{j},RU{j}),min(LD_min{j},RD{j}),min(LU_min{j},LD{j}));
 [RU_min3{j},RD_min3{j},LD_min3{j},LU_min3{j}] = deal(min(RU_min{j},RD{j}),min(RD_min{j},LD{j}),min(LD_min{j},LU{j}),min(LU_min{j},RU{j}));
 [RU_max2{j},RD_max2{j},LD_max2{j},LU_max2{j}] = deal(max(RU_max{j},LU{j}),max(RD_max{j},RU{j}),max(LD_max{j},RD{j}),max(LU_max{j},LD{j}));
 [RU_max3{j},RD_max3{j},LD_max3{j},LU_max3{j}] = deal(max(RU_max{j},RD{j}),max(RD_max{j},LD{j}),max(LD_max{j},LU{j}),max(LU_max{j},RU{j}));
end,    clear RU_min RD_min LU_min LD_min RU_max RD_max LU_max LD_max
 g.min48c = Cooc1split(vercat(vercat(RU_min2,LD_min2),vercat(RD_min3,LU_min3)),3,'col',T) +...
     Cooc1split(horcat(horcat(RD_min2,LU_min2),horcat(RU_min3,LD_min3)),3,'col',T);
 g.max48c = Cooc1split(vercat(vercat(RU_max2,LD_max2),vercat(RD_max3,LU_max3)),3,'col',T) +...
     Cooc1split(horcat(horcat(RD_max2,LU_max2),horcat(RU_max3,LD_max3)),3,'col',T);

 % minmax54 -- to be symmetrized as mnmx, directional, hv-symmetrical
for j=1:three
 [RU_min4{j},RD_min4{j},LD_min4{j},LU_min4{j}] = deal(min(RU_min2{j},RD{j}),min(RD_min2{j},LD{j}),min(LD_min2{j},LU{j}),min(LU_min2{j},RU{j}));
 [RU_min5{j},RD_min5{j},LD_min5{j},LU_min5{j}] = deal(min(RU_min3{j},LU{j}),min(RD_min3{j},RU{j}),min(LD_min3{j},RD{j}),min(LU_min3{j},LD{j}));
 [RU_max4{j},RD_max4{j},LD_max4{j},LU_max4{j}] = deal(max(RU_max2{j},RD{j}),max(RD_max2{j},LD{j}),max(LD_max2{j},LU{j}),max(LU_max2{j},RU{j}));
 [RU_max5{j},RD_max5{j},LD_max5{j},LU_max5{j}] = deal(max(RU_max3{j},LU{j}),max(RD_max3{j},RU{j}),max(LD_max3{j},RD{j}),max(LU_max3{j},LD{j}));
end,    clear RU_min2 RD_min2 LD_min2 LU_min2 RU_min3 RD_min3 LD_min3 LU_min3 
        clear RU_max2 RD_max2 LD_max2 LU_max2 RU_max3 RD_max3 LD_max3 LU_max3 
 g.min54c = Cooc1split(vercat(vercat(RU_min4,LD_min4),vercat(RD_min5,LU_min5)),3,'col',T) +...
     Cooc1split(horcat(horcat(RD_min4,LU_min4),horcat(RU_min5,LD_min5)),3,'col',T);
 g.max54c = Cooc1split(vercat(vercat(RU_max4,LD_max4),vercat(RD_max5,LU_max5)),3,'col',T) +...
     Cooc1split(horcat(horcat(RD_max4,LU_max4),horcat(RU_max5,LD_max5)),3,'col',T);
 
%**** FUNCTION ****%
function C = vercat(A,B)
if iscell(A) & iscell(B)
    C{1} = [A{1};B{1}];
    for j=2:numel(A)
        C{j} = [A{j};B{j}];
    end
end

%**** FUNCTION ****%
function C = horcat(A,B)
if iscell(A) & iscell(B)
    C{1} = [A{1},B{1}];
    for j=2:numel(A)
        C{j} = [A{j},B{j}];
    end
end

%**** FUNCTION ****%
function g = all2nd(X,q)
%
% X must be a matrix of doubles or singles (the image) and q is the 
% quantization step (any positive number).
%
% Recommended values of q are c, 1.5c, 2c, where c is the central
% coefficient in the differential (at X(I,J)).
%
% This function outputs co-occurrences of ALL 2nd-order residuals
% listed in Figure 1 in our journal HUGO paper (version from June 14), 
% including the naming convention.
%
% List of outputted features:
%
% 1a) spam12h
% 1b) spam12v (orthogonal-spam)
% 1c) minmax21
% 1d) minmax41
% 1e) minmax24h (24v is also outputted but not listed in Figure 1)
% 1f) minmax32
%
% Naming convention:
%
% name = {type}{f}{sigma}{scan}
% type \in {spam, minmax}
% f \in {1,2,3,4,5} number of filters that are "minmaxed"
% sigma \in {1,2,3,4,8} symmetry index
% scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
% h and v scans).
%
% All even residuals are implemented the same way simply by
% narrowing the range for I and J and replacing the residuals.
%
% Note1: The term X(I,J) should always have the "-" sign.
% Note2: This script does not include s, so, cout, cin versions (weak).
%
% This function calls Residual.m, Cooc1split.m, and Quant.m

[M N three] = size(X); [T,order] = deal(2,4);
% 2nd-order residuals are implemented using Residual.m
for j=1:three
 [Dh{j},Dv{j},Dd{j},Dm{j}] = deal(Residual(X(:,:,j),2,'hor'),Residual(X(:,:,j),2,'ver'),Residual(X(:,:,j),2,'diag'),Residual(X(:,:,j),2,'mdiag'));
 [Yh{j},Yv{j},Yd{j},Ym{j}] = deal(Quant(Dh{j},q,T),Quant(Dv{j},q,T),Quant(Dd{j},q,T),Quant(Dm{j},q,T));
end,    clear Dh Dv Dd Dm X

% spam12h/v
 g.spam12h = Cooc1split(Yh,order,'hor',T) + Cooc1split(Yv,order,'ver',T);
 g.spam12v = Cooc1split(Yh,order,'ver',T) + Cooc1split(Yv,order,'hor',T);

% minmax21
for j=1:three
 [Dmin{j},Dmax{j}] = deal(min(Yh{j},Yv{j}),max(Yh{j},Yv{j}));
end
 g.min21 = Cooc1split(Dmin,order,'hor',T) + Cooc1split(Dmin,order,'ver',T);
 g.max21 = Cooc1split(Dmax,order,'hor',T) + Cooc1split(Dmax,order,'ver',T);

% minmax41   
for j=1:three
 [Dmin2{j},Dmax2{j}] = deal(min(Dmin{j},min(Yd{j},Ym{j})),max(Dmax{j},max(Yd{j},Ym{j})));
end
 g.min41 = Cooc1split(Dmin2,order,'hor',T) + Cooc1split(Dmin2,order,'ver',T);
 g.max41 = Cooc1split(Dmax2,order,'hor',T) + Cooc1split(Dmax2,order,'ver',T);
   	clear Dmin2 Dmax2

 % minmax24h,v -- both "not bad," h slightly better, directional, hv-nonsymmetrical, to be symmetrized as mnmx
for j=1:three
 [RU_min2{j},RD_min2{j},RU_min3{j},LU_min3{j}] = deal(min(Ym{j},Yh{j}),min(Yd{j},Yh{j}),min(Ym{j},Yv{j}),min(Yd{j},Yv{j}));
 [RU_max2{j},RD_max2{j},RU_max3{j},LU_max3{j}] = deal(max(Ym{j},Yh{j}),max(Yd{j},Yh{j}),max(Ym{j},Yv{j}),max(Yd{j},Yv{j}));
end,    clear Yh Yv
 g.min24h = Cooc1split(vercat(RU_min2,RD_min2),order,'hor',T) + Cooc1split(horcat(RU_min3,LU_min3),order,'ver',T);
 g.min24v = Cooc1split(horcat(RU_min2,RD_min2),order,'ver',T) + Cooc1split(vercat(RU_min3,LU_min3),order,'hor',T);
 g.max24h = Cooc1split(vercat(RU_max2,RD_max2),order,'hor',T) + Cooc1split(horcat(RU_max3,LU_max3),order,'ver',T);
 g.max24v = Cooc1split(horcat(RU_max2,RD_max2),order,'ver',T) + Cooc1split(vercat(RU_max3,LU_max3),order,'hor',T);

 % minmax32 -- good, directional, hv-symmetrical, to be symmetrized as mnmx
for j=1:three
 [RU_min{j},RD_min{j}] = deal(min(Dmin{j},Ym{j}),min(Dmin{j},Yd{j})); 
 [RU_max{j},RD_max{j}] = deal(max(Dmax{j},Ym{j}),max(Dmax{j},Yd{j})); 
end,    clear Dmin Dmax
 g.min32 = Cooc1split(vercat(RU_min,RD_min),order,'hor',T) + Cooc1split(horcat(RU_min,RD_min),order,'ver',T);
 g.max32 = Cooc1split(vercat(RU_max,RD_max),order,'hor',T) + Cooc1split(horcat(RU_max,RD_max),order,'ver',T);

%**** FUNCTION ****%
function g = all2ndc(X,q,T)
%
% X must be a matrix of doubles or singles (the image) and q is the 
% quantization step (any positive number).
%
% Recommended values of q are c, 1.5c, 2c, where c is the central
% coefficient in the differential (at X(I,J)).
%
% This function outputs co-occurrences of ALL 2nd-order residuals
% listed in Figure 1 in our journal HUGO paper (version from June 14), 
% including the naming convention.
%
% List of outputted features:
%
% 1a) spam12h
% 1b) spam12v (orthogonal-spam)
% 1c) minmax21
% 1d) minmax41
% 1e) minmax24h (24v is also outputted but not listed in Figure 1)
% 1f) minmax32
%
% Naming convention:
%
% name = {type}{f}{sigma}{scan}
% type \in {spam, minmax}
% f \in {1,2,3,4,5} number of filters that are "minmaxed"
% sigma \in {1,2,3,4,8} symmetry index
% scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
% h and v scans).
%
% All even residuals are implemented the same way simply by
% narrowing the range for I and J and replacing the residuals.
%
% Note1: The term X(I,J) should always have the "-" sign.
% Note2: This script does not include s, so, cout, cin versions (weak).
%
% This function calls Residual.m, Cooc1split.m, and Quant.m

[M N three] = size(X); order = 4;
% 2nd-order residuals are implemented using Residual.m
for j=1:three
 [Dh{j},Dv{j},Dd{j},Dm{j}] = deal(Residual(X(:,:,j),2,'hor'),Residual(X(:,:,j),2,'ver'),Residual(X(:,:,j),2,'diag'),Residual(X(:,:,j),2,'mdiag'));
 [Yh{j},Yv{j},Yd{j},Ym{j}] = deal(Quant(Dh{j},q,T),Quant(Dv{j},q,T),Quant(Dd{j},q,T),Quant(Dm{j},q,T));
end,    clear Dh Dv Dd Dm X
% spam12h/v
 g.spam12c = Cooc1split(Yh,3,'col',T) + Cooc1split(Yv,3,'col',T);

% minmax21
for j=1:three
 [Dmin{j},Dmax{j}] = deal(min(Yh{j},Yv{j}),max(Yh{j},Yv{j}));
end
 g.min21c = Cooc1split(Dmin,3,'col',T);
 g.max21c = Cooc1split(Dmax,3,'col',T);

% minmax41   
for j=1:three
 [Dmin2{j},Dmax2{j}] = deal(min(Dmin{j},min(Yd{j},Ym{j})),max(Dmax{j},max(Yd{j},Ym{j})));
end
 g.min41c = Cooc1split(Dmin2,3,'col',T);      clear Dmin2
 g.max41c = Cooc1split(Dmax2,3,'col',T);      clear Dmax2

 % minmax24h,v -- both "not bad," h slightly better, directional, hv-nonsymmetrical, to be symmetrized as mnmx
for j=1:three
 [RU_min2{j},RD_min2{j},RU_min3{j},LU_min3{j}] = deal(min(Ym{j},Yh{j}),min(Yd{j},Yh{j}),min(Ym{j},Yv{j}),min(Yd{j},Yv{j}));
 [RU_max2{j},RD_max2{j},RU_max3{j},LU_max3{j}] = deal(max(Ym{j},Yh{j}),max(Yd{j},Yh{j}),max(Ym{j},Yv{j}),max(Yd{j},Yv{j}));
end,    clear Yh Yv
 g.min24c = Cooc1split(vercat(RU_min2,RD_min2),3,'col',T) + Cooc1split(horcat(RU_min3,LU_min3),3,'col',T);
 g.max24c = Cooc1split(vercat(RU_max2,RD_max2),3,'col',T) + Cooc1split(horcat(RU_max3,LU_max3),3,'col',T);

 % minmax32 -- good, directional, hv-symmetrical, to be symmetrized as mnmx
for j=1:three
 [RU_min{j},RD_min{j}] = deal(min(Dmin{j},Ym{j}),min(Dmin{j},Yd{j})); 
 [RU_max{j},RD_max{j}] = deal(max(Dmax{j},Ym{j}),max(Dmax{j},Yd{j})); 
end,    clear Dmin Dmax
 g.min32c = Cooc1split(vercat(RU_min,RD_min),3,'col',T);
 g.max32c = Cooc1split(vercat(RU_max,RD_max),3,'col',T);

%**** FUNCTION ****%
function g = all3x3(X,q)
% This function outputs co-occurrences of ALL residuals based on the
% KB kernel and its "halves" (EDGE residuals) as listed in Figure 1
% in our journal HUGO paper (version from June 14), including the naming
% convention.

[M N three] = size(X);   [T,order] = deal(2,4); 
% spam11 (old name KB residual), good, non-directional, hv-symmetrical, to be symmetrized as spam
for j=1:three
 D{j} = Residual(X(:,:,j),2,'KB'); Y{j} = Quant(D{j},q,T);
end
 g.spam11 = Cooc1split(Y,order,'hor',T) + Cooc1split(Y,order,'ver',T);

% EDGE residuals
for j=1:three
 D{j} = Residual(X(:,:,j),2,'edge-h'); Du{j} = D{j}(:,1:size(D{j},2)/2); Db{j} = D{j}(:,size(D{j},2)/2+1:end);
 D{j} = Residual(X(:,:,j),2,'edge-v'); Dl{j} = D{j}(:,1:size(D{j},2)/2); Dr{j} = D{j}(:,size(D{j},2)/2+1:end);
 [Yu{j},Yb{j},Yl{j},Yr{j}] = deal(Quant(Du{j},q,T),Quant(Db{j},q,T),Quant(Dl{j},q,T),Quant(Dr{j},q,T));
end,     clear D

% spam14h,v  not bad, directional, hv-nonsym, to be symmetrized as spam
 g.spam14v = Cooc1split(horcat(Yu,Yb),order,'ver',T) + Cooc1split(vercat(Yl,Yr),order,'hor',T);
 g.spam14h = Cooc1split(vercat(Yu,Yb),order,'hor',T) + Cooc1split(horcat(Yl,Yr),order,'ver',T);

% minmax24 -- EXCELLENT, directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin1{j},Dmin2{j},Dmin3{j},Dmin4{j}] = deal(min(Yu{j},Yl{j}),min(Yb{j},Yr{j}),min(Yu{j},Yr{j}),min(Yb{j},Yl{j}));
 [Dmax1{j},Dmax2{j},Dmax3{j},Dmax4{j}] = deal(max(Yu{j},Yl{j}),max(Yb{j},Yr{j}),max(Yu{j},Yr{j}),max(Yb{j},Yl{j}));
end 
 g.min24 = Cooc1split(horcat(horcat(Dmin1,Dmin2),horcat(Dmin3,Dmin4)),order,'ver',T) + ...
     Cooc1split(vercat(vercat(Dmin1,Dmin2),vercat(Dmin3,Dmin4)),order,'hor',T);
 g.max24 = Cooc1split(horcat(horcat(Dmax1,Dmax2),horcat(Dmax3,Dmax4)),order,'ver',T) + ...
     Cooc1split(vercat(vercat(Dmax1,Dmax2),vercat(Dmax3,Dmax4)),order,'hor',T);
    clear Dmin3 Dmin4 Dmax3 Dmax4
    
% minmax22 - hv-nonsymmetrical
% min22h -- good, to be symmetrized as mnmx, directional, hv-nonsymmetrical
% min22v -- EXCELLENT - to be symmetrized as mnmx, directional,
for j=1:three
 [UEq_min{j},REq_min{j}] = deal(min(Yu{j},Yb{j}),min(Yr{j},Yl{j}));
 [UEq_max{j},REq_max{j}] = deal(max(Yu{j},Yb{j}),max(Yr{j},Yl{j}));
end
 g.min22h = Cooc1split(UEq_min,order,'hor',T) + Cooc1split(REq_min,order,'ver',T);
 g.min22v = Cooc1split(UEq_min,order,'ver',T) + Cooc1split(REq_min,order,'hor',T);
 g.max22h = Cooc1split(UEq_max,order,'hor',T) + Cooc1split(REq_max,order,'ver',T);
 g.max22v = Cooc1split(UEq_max,order,'ver',T) + Cooc1split(REq_max,order,'hor',T);

% minmax41 -- good, non-directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin5{j},Dmax5{j}] = deal(min(Dmin1{j},Dmin2{j}),max(Dmax1{j},Dmax2{j}));
end,     clear Dmin1 Dmin2 Dmax1 Dmax2 
 g.min41 = Cooc1split(Dmin5,order,'ver',T) + Cooc1split(Dmin5,order,'hor',T);
 g.max41 = Cooc1split(Dmax5,order,'ver',T) + Cooc1split(Dmax5,order,'hor',T);

%**** FUNCTION ****%
function g = all3x3c(X,q,T)
% This function outputs co-occurrences of ALL residuals based on the
% KB kernel and its "halves" (EDGE residuals) as listed in Figure 1
% in our journal HUGO paper (version from June 14), including the naming
% convention.

[M N three] = size(X);   order = 4;
% spam11 (old name KB residual), good, non-directional, hv-symmetrical, to be symmetrized as spam
for j=1:three
 D{j} = Residual(X(:,:,j),2,'KB'); Y{j} = Quant(D{j},q,T);
end
 g.spam11c = Cooc1split(Y,3,'col',T);

% EDGE residuals
for j=1:three
 D{j} = Residual(X(:,:,j),2,'edge-h'); Du{j} = D{j}(:,1:size(D{j},2)/2); Db{j} = D{j}(:,size(D{j},2)/2+1:end);
 D{j} = Residual(X(:,:,j),2,'edge-v'); Dl{j} = D{j}(:,1:size(D{j},2)/2); Dr{j} = D{j}(:,size(D{j},2)/2+1:end);
 [Yu{j},Yb{j},Yl{j},Yr{j}] = deal(Quant(Du{j},q,T),Quant(Db{j},q,T),Quant(Dl{j},q,T),Quant(Dr{j},q,T));
end,     clear D

% spam14h,v  not bad, directional, hv-nonsym, to be symmetrized as spam
 g.spam14c = Cooc1split(vercat(Yu,Yb),3,'col',T) + Cooc1split(horcat(Yl,Yr),3,'col',T);

% minmax24 -- EXCELLENT, directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin1{j},Dmin2{j},Dmin3{j},Dmin4{j}] = deal(min(Yu{j},Yl{j}),min(Yb{j},Yr{j}),min(Yu{j},Yr{j}),min(Yb{j},Yl{j}));
 [Dmax1{j},Dmax2{j},Dmax3{j},Dmax4{j}] = deal(max(Yu{j},Yl{j}),max(Yb{j},Yr{j}),max(Yu{j},Yr{j}),max(Yb{j},Yl{j}));
end 
 g.min24c = Cooc1split(vercat(vercat(Dmin1,Dmin2),vercat(Dmin3,Dmin4)),3,'col',T);
 g.max24c = Cooc1split(vercat(vercat(Dmax1,Dmax2),vercat(Dmax3,Dmax4)),3,'col',T);
    clear Dmin3 Dmin4 Dmax3 Dmax4
    
% minmax22 - hv-nonsymmetrical
% min22h -- good, to be symmetrized as mnmx, directional, hv-nonsymmetrical
% min22v -- EXCELLENT - to be symmetrized as mnmx, directional,
for j=1:three
 [UEq_min{j},REq_min{j}] = deal(min(Yu{j},Yb{j}),min(Yr{j},Yl{j}));
 [UEq_max{j},REq_max{j}] = deal(max(Yu{j},Yb{j}),max(Yr{j},Yl{j}));
end
 g.min22c = Cooc1split(UEq_min,3,'col',T) + Cooc1split(REq_min,3,'col',T);
 g.max22c = Cooc1split(UEq_max,3,'col',T) + Cooc1split(REq_max,3,'col',T);

% minmax41 -- good, non-directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin5{j},Dmax5{j}] = deal(min(Dmin1{j},Dmin2{j}),max(Dmax1{j},Dmax2{j}));
end,     clear Dmin1 Dmin2 Dmax1 Dmax2 
 g.min41c = Cooc1split(Dmin5,3,'col',T);
 g.max41c = Cooc1split(Dmax5,3,'col',T);

%**** FUNCTION ****%
function g = all5x5(X,q)
% This function outputs co-occurrences of ALL residuals based on the
% KV kernel and its "halves" (EDGE residuals) as listed in Figure 1
% in our journal HUGO paper (version from June 14), including the naming
% convention.
[M N three] = size(X);  [I,J,T,order] = deal(3:M-2,3:N-2,2,4);

% spam11 (old name KV residual), good, non-directional, hv-symmetrical, to be symmetrized as spam
for j=1:three
 D{j} = Residual(X(:,:,j),3,'KV'); Y{j} = Quant(D{j},q,T);
end
 g.spam11 = Cooc1split(Y,order,'hor',T) + Cooc1split(Y,order,'ver',T);

% EDGE residuals    
for j=1:three
 Du = 8*X(I,J-1,j)+8*X(I-1,J,j)+8*X(I,J+1,j)-6*X(I-1,J-1,j)-6*X(I-1,J+1,j)-2*X(I,J-2,j)-2*X(I,J+2,j)-2*X(I-2,J,j)+2*X(I-1,J-2,j)+2*X(I-2,J-1,j)+2*X(I-2,J+1,j)+2*X(I-1,J+2,j)-X(I-2,J-2,j)-X(I-2,J+2,j)-12*X(I,J,j);
 Dr = 8*X(I-1,J,j)+8*X(I,J+1,j)+8*X(I+1,J,j)-6*X(I-1,J+1,j)-6*X(I+1,J+1,j)-2*X(I-2,J,j)-2*X(I+2,J,j)-2*X(I,J+2,j)+2*X(I-2,J+1,j)+2*X(I-1,J+2,j)+2*X(I+1,J+2,j)+2*X(I+2,J+1,j)-X(I-2,J+2,j)-X(I+2,J+2,j)-12*X(I,J,j);
 Db = 8*X(I,J+1,j)+8*X(I+1,J,j)+8*X(I,J-1,j)-6*X(I+1,J+1,j)-6*X(I+1,J-1,j)-2*X(I,J-2,j)-2*X(I,J+2,j)-2*X(I+2,J,j)+2*X(I+1,J+2,j)+2*X(I+2,J+1,j)+2*X(I+2,J-1,j)+2*X(I+1,J-2,j)-X(I+2,J+2,j)-X(I+2,J-2,j)-12*X(I,J,j);
 Dl = 8*X(I+1,J,j)+8*X(I,J-1,j)+8*X(I-1,J,j)-6*X(I+1,J-1,j)-6*X(I-1,J-1,j)-2*X(I-2,J,j)-2*X(I+2,J,j)-2*X(I,J-2,j)+2*X(I+2,J-1,j)+2*X(I+1,J-2,j)+2*X(I-1,J-2,j)+2*X(I-2,J-1,j)-X(I+2,J-2,j)-X(I-2,J-2,j)-12*X(I,J,j);
 [Yu{j},Yb{j},Yl{j},Yr{j}] = deal(Quant(Du,q,T),Quant(Db,q,T),Quant(Dl,q,T),Quant(Dr,q,T));
end,    clear Du Db Dl Dr

% spam14v  not bad, directional, hv-nonsym, to be symmetrized as spam
 g.spam14v = Cooc1split(horcat(Yu,Yb),order,'ver',T) + Cooc1split(vercat(Yl,Yr),order,'hor',T);
 g.spam14h = Cooc1split(vercat(Yu,Yb),order,'hor',T) + Cooc1split(horcat(Yl,Yr),order,'ver',T);

% minmax24 -- EXCELLENT, directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin1{j},Dmin2{j},Dmin3{j},Dmin4{j}] = deal(min(Yu{j},Yl{j}),min(Yb{j},Yr{j}),min(Yu{j},Yr{j}),min(Yb{j},Yl{j}));
 [Dmax1{j},Dmax2{j},Dmax3{j},Dmax4{j}] = deal(max(Yu{j},Yl{j}),max(Yb{j},Yr{j}),max(Yu{j},Yr{j}),max(Yb{j},Yl{j}));
end
 g.min24 = Cooc1split(horcat(horcat(Dmin1,Dmin2),horcat(Dmin3,Dmin4)),order,'ver',T) + ...
     Cooc1split(vercat(vercat(Dmin1,Dmin2),vercat(Dmin3,Dmin4)),order,'hor',T);
 g.max24 = Cooc1split(horcat(horcat(Dmax1,Dmax2),horcat(Dmax3,Dmax4)),order,'ver',T) + ...
     Cooc1split(vercat(vercat(Dmax1,Dmax2),vercat(Dmax3,Dmax4)),order,'hor',T);
    clear Dmin3 Dmax3 Dmin4 Dmax4
% minmax22 - hv-nonsymmetrical
% min22h -- good, to be symmetrized as mnmx, directional, hv-nonsymmetrical
% min22v -- EXCELLENT - to be symmetrized as mnmx, directional,
for j=1:three
 [UEq_min{j},REq_min{j}] = deal(min(Yu{j},Yb{j}),min(Yr{j},Yl{j}));
 [UEq_max{j},REq_max{j}] = deal(max(Yu{j},Yb{j}),max(Yr{j},Yl{j}));
end
 g.min22h = Cooc1split(UEq_min,order,'hor',T) + Cooc1split(REq_min,order,'ver',T);
 g.min22v = Cooc1split(UEq_min,order,'ver',T) + Cooc1split(REq_min,order,'hor',T);
 g.max22h = Cooc1split(UEq_max,order,'hor',T) + Cooc1split(REq_max,order,'ver',T);
 g.max22v = Cooc1split(UEq_max,order,'ver',T) + Cooc1split(REq_max,order,'hor',T);

% minmax41 -- good, non-directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin5{j},Dmax5{j}] = deal(min(Dmin1{j},Dmin2{j}),max(Dmax1{j},Dmax2{j}));
end
 g.min41 = Cooc1split(Dmin5,order,'ver',T) + Cooc1split(Dmin5,order,'hor',T);
 g.max41 = Cooc1split(Dmax5,order,'ver',T) + Cooc1split(Dmax5,order,'hor',T);

%**** FUNCTION ****%
function g = all5x5c(X,q,T)
% This function outputs co-occurrences of ALL residuals based on the
% KV kernel and its "halves" (EDGE residuals) as listed in Figure 1
% in our journal HUGO paper (version from June 14), including the naming
% convention.
[M N three] = size(X);  [I,J,order] = deal(3:M-2,3:N-2,4);
% spam11 (old name KV residual), good, non-directional, hv-symmetrical, to be symmetrized as spam
for j=1:three
 D{j} = Residual(X(:,:,j),3,'KV'); Y{j} = Quant(D{j},q,T);
end
 g.spam11c = Cooc1split(Y,3,'col',T);

% EDGE residuals    
for j=1:three
 Du = 8*X(I,J-1,j)+8*X(I-1,J,j)+8*X(I,J+1,j)-6*X(I-1,J-1,j)-6*X(I-1,J+1,j)-2*X(I,J-2,j)-2*X(I,J+2,j)-2*X(I-2,J,j)+2*X(I-1,J-2,j)+2*X(I-2,J-1,j)+2*X(I-2,J+1,j)+2*X(I-1,J+2,j)-X(I-2,J-2,j)-X(I-2,J+2,j)-12*X(I,J,j);
 Dr = 8*X(I-1,J,j)+8*X(I,J+1,j)+8*X(I+1,J,j)-6*X(I-1,J+1,j)-6*X(I+1,J+1,j)-2*X(I-2,J,j)-2*X(I+2,J,j)-2*X(I,J+2,j)+2*X(I-2,J+1,j)+2*X(I-1,J+2,j)+2*X(I+1,J+2,j)+2*X(I+2,J+1,j)-X(I-2,J+2,j)-X(I+2,J+2,j)-12*X(I,J,j);
 Db = 8*X(I,J+1,j)+8*X(I+1,J,j)+8*X(I,J-1,j)-6*X(I+1,J+1,j)-6*X(I+1,J-1,j)-2*X(I,J-2,j)-2*X(I,J+2,j)-2*X(I+2,J,j)+2*X(I+1,J+2,j)+2*X(I+2,J+1,j)+2*X(I+2,J-1,j)+2*X(I+1,J-2,j)-X(I+2,J+2,j)-X(I+2,J-2,j)-12*X(I,J,j);
 Dl = 8*X(I+1,J,j)+8*X(I,J-1,j)+8*X(I-1,J,j)-6*X(I+1,J-1,j)-6*X(I-1,J-1,j)-2*X(I-2,J,j)-2*X(I+2,J,j)-2*X(I,J-2,j)+2*X(I+2,J-1,j)+2*X(I+1,J-2,j)+2*X(I-1,J-2,j)+2*X(I-2,J-1,j)-X(I+2,J-2,j)-X(I-2,J-2,j)-12*X(I,J,j);
 [Yu{j},Yb{j},Yl{j},Yr{j}] = deal(Quant(Du,q,T),Quant(Db,q,T),Quant(Dl,q,T),Quant(Dr,q,T));
end,    clear Du Db Dl Dr

% spam14v  not bad, directional, hv-nonsym, to be symmetrized as spam
 g.spam14c = Cooc1split(horcat(Yu,Yb),3,'col',T) + Cooc1split(vercat(Yl,Yr),3,'col',T);

% minmax24 -- EXCELLENT, directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin1{j},Dmin2{j},Dmin3{j},Dmin4{j}] = deal(min(Yu{j},Yl{j}),min(Yb{j},Yr{j}),min(Yu{j},Yr{j}),min(Yb{j},Yl{j}));
 [Dmax1{j},Dmax2{j},Dmax3{j},Dmax4{j}] = deal(max(Yu{j},Yl{j}),max(Yb{j},Yr{j}),max(Yu{j},Yr{j}),max(Yb{j},Yl{j}));
end
 g.min24c = Cooc1split(horcat(horcat(Dmin1,Dmin2),horcat(Dmin3,Dmin4)),3,'col',T);      clear Dmin3 Dmin4 
 g.max24c = Cooc1split(horcat(horcat(Dmax1,Dmax2),horcat(Dmax3,Dmax4)),3,'col',T);      clear Dmax3 Dmax4
 
% minmax22 - hv-nonsymmetrical
% min22h -- good, to be symmetrized as mnmx, directional, hv-nonsymmetrical
% min22v -- EXCELLENT - to be symmetrized as mnmx, directional,
for j=1:three
 [UEq_min{j},REq_min{j}] = deal(min(Yu{j},Yb{j}),min(Yr{j},Yl{j}));
 [UEq_max{j},REq_max{j}] = deal(max(Yu{j},Yb{j}),max(Yr{j},Yl{j}));
end
 g.min22c = Cooc1split(UEq_min,3,'col',T) + Cooc1split(REq_min,3,'col',T);        clear UEq_min REq_min
 g.max22c = Cooc1split(UEq_max,3,'col',T) + Cooc1split(REq_max,3,'col',T);        clear UEq_max REq_max

% minmax41 -- good, non-directional, hv-sym, to be symmetrized as mnmx
for j=1:three
 [Dmin5{j},Dmax5{j}] = deal(min(Dmin1{j},Dmin2{j}),max(Dmax1{j},Dmax2{j}));
end
 g.min41c = Cooc1split(Dmin5,3,'col',T);
 g.max41c = Cooc1split(Dmax5,3,'col',T);
 
 
%**** FUNCTION ****%
function f = Cooc1split(D,order,type,T)
% Concatenate two sets of co-occurrence features from the residual split to two parts 

for j=1:3
    D_green{j}   = [D{j}(2:2:end,1:2:end),D{j}(1:2:end,2:2:end)];
    D_redblue{j} = [D{j}(1:2:end,1:2:end),D{j}(2:2:end,2:2:end)];
end

switch type
case 'hor'
    for j=1:3
        fj{j} = [Cooc1(D_green{j},order,type,T); Cooc1(D_redblue{j},order,type,T)];
    end
    f = fj{1}+fj{2}+fj{3};
case 'ver'
    for j=1:3
        fj{j} = [Cooc1(D_green{j},order,type,T); Cooc1(D_redblue{j},order,type,T)];
    end
    f = fj{1}+fj{2}+fj{3};
case 'col'
    f = [Cooc1(D_green,order,type,T); Cooc1(D_redblue,order,type,T)];
end

%**** FUNCTION ****%
function f = Cooc1(D,order,type,T)
% Co-occurrence operator to be appied to a 2D array of residuals D \in [-T,T]
% T     ... threshold
% order ... cooc order \in {1,2,3,4,5}
% type  ... cooc type \in {hor,ver,diag,mdiag,square,square-ori,hvdm}
% f     ... an array of length (2T+1)^order  (1D array!)

if iscell(D),
  for j=1:numel(D)
    if max(abs(D{j}(:))) > T, fprintf('*** ERROR in Cooc1.m: Residual out of range ***\n'), end
  end  
else
    if max(abs(D(:))) > T, fprintf('*** ERROR in Cooc1.m: Residual out of range ***\n'), end
end
B = 2*T+1;  z = (B^order-1)/2;  range = (-z:z);
switch order
    case 1
        f = histc(D(:),-T:T);
    case 2
        f = zeros(B,B);
        if strcmp(type,'hor'),   L = D(:,1:end-1); R = D(:,2:end);end
        if strcmp(type,'ver'),   L = D(1:end-1,:); R = D(2:end,:);end
        if strcmp(type,'diag'),  L = D(1:end-1,1:end-1); R = D(2:end,2:end);end
        if strcmp(type,'mdiag'), L = D(1:end-1,2:end); R = D(2:end,1:end-1);end
        A = L(:)+B*R(:);   clear L R
        f = histc(A,range); 
    case 3
        f = zeros(B,B,B);
        if strcmp(type,'hor'),   L = D(:,1:end-2); C = D(:,2:end-1); R = D(:,3:end);end
        if strcmp(type,'ver'),   L = D(1:end-2,:); C = D(2:end-1,:); R = D(3:end,:);end
        if strcmp(type,'diag'),  L = D(1:end-2,1:end-2); C = D(2:end-1,2:end-1); R = D(3:end,3:end);end
        if strcmp(type,'mdiag'), L = D(1:end-2,3:end); C = D(2:end-1,2:end-1); R = D(3:end,1:end-2);end
        if strcmp(type,'col'),   L = D{1}; C = D{2}; R = D{3};  end       
        A = L(:)+B*C(:)+B^2*R(:);    clear L C R
        f = histc(A,range); 
    case 4
        f = zeros(B,B,B,B);
        if strcmp(type,'hor'),    L = D(:,1:end-3); C = D(:,2:end-2); E = D(:,3:end-1); R = D(:,4:end);end
        if strcmp(type,'ver'),    L = D(1:end-3,:); C = D(2:end-2,:); E = D(3:end-1,:); R = D(4:end,:);end
        if strcmp(type,'diag'),   L = D(1:end-3,1:end-3); C = D(2:end-2,2:end-2); E = D(3:end-1,3:end-1); R = D(4:end,4:end);end
        if strcmp(type,'mdiag'),  L = D(4:end,1:end-3); C = D(3:end-1,2:end-2); E = D(2:end-2,3:end-1); R = D(1:end-3,4:end);end
        if strcmp(type,'square'), L = D(2:end,1:end-1); C = D(2:end,2:end); E = D(1:end-1,2:end); R = D(1:end-1,1:end-1);end
        if strcmp(type,'square-ori'), [M, N] = size(D); Dh = D(:,1:M); Dv = D(:,M+1:2*M);
                                  L = Dh(2:end,1:end-1); C = Dv(2:end,2:end); E = Dh(1:end-1,2:end); R = Dv(1:end-1,1:end-1);end
        if strcmp(type,'hvdm'),   [M, N] = size(D); L = D(:,1:M); C = D(:,M+1:2*M); E = D(:,2*M+1:3*M); R = D(:,3*M+1:4*M);end
        if strcmp(type,'s-in'),   [M, N] = size(D); Dh = D(:,1:M); Dv = D(:,M+1:2*M); Dh1 = D(:,2*M+1:3*M); Dv1 = D(:,3*M+1:4*M);
                                  L = Dh(2:end,1:end-1); C = Dh1(2:end,2:end); E = Dh1(1:end-1,2:end); R = Dh(1:end-1,1:end-1);end
        if strcmp(type,'s-out'),  [M, N] = size(D); Dh = D(:,1:M); Dv = D(:,M+1:2*M); Dh1 = D(:,2*M+1:3*M); Dv1 = D(:,3*M+1:4*M);
                                  L = Dh1(2:end,1:end-1); C = Dh(2:end,2:end); E = Dh(1:end-1,2:end); R = Dh1(1:end-1,1:end-1);end
        if strcmp(type,'ori-in'), [M, N] = size(D); Dh = D(:,1:M); Dv = D(:,M+1:2*M); Dh1 = D(:,2*M+1:3*M); Dv1 = D(:,3*M+1:4*M);
                                  L = Dh(2:end,1:end-1); C = Dv1(2:end,2:end); E = Dh1(1:end-1,2:end); R = Dv(1:end-1,1:end-1);end
        if strcmp(type,'ori-out'),[M, N] = size(D); Dh = D(:,1:M); Dv = D(:,M+1:2*M); Dh1 = D(:,2*M+1:3*M); Dv1 = D(:,3*M+1:4*M);
                                  L = Dh1(2:end,1:end-1); C = Dv(2:end,2:end); E = Dh(1:end-1,2:end); R = Dv1(1:end-1,1:end-1);end
        A = L(:)+B*C(:)+B^2*E(:)+B^3*R(:);    clear L C E R
        f = histc(A,range); 
    case 5
        f = zeros(B,B,B,B,B);
        if strcmp(type,'hor'),L = D(:,1:end-4); C = D(:,2:end-3); E = D(:,3:end-2); F = D(:,4:end-1); R = D(:,5:end);end
        if strcmp(type,'ver'),L = D(1:end-4,:); C = D(2:end-3,:); E = D(3:end-2,:); F = D(4:end-1,:); R = D(5:end,:);end
        A = L(:)+B*C(:)+B^2*E(:)+B^3*F(:)+B^4*R(:);    clear L C E F R
        f = histc(A,range); 
end

% normalization
f = double(f);
fsum = sum(f(:));
if fsum>0, f = f/fsum; end

%**** FUNCTION ****%
function A4 = CFAsplit(A,dim)
% A     ... n by 1 array containing some quantity computed from sliding blocks in matrix L 
% A4    ... 4 cells containing distinct subsets of A. In one subset,
% sliding blocks have the same pixel type (referring to Bayer CFA) in their upper left corner

block_numbers = 1:prod(dim);    % enumerate the blocks
block_numbers = reshape(block_numbers,dim);

type{1} = block_numbers(1:2:end,1:2:end);
type{2} = block_numbers(2:2:end,1:2:end);
type{3} = block_numbers(1:2:end,2:2:end);
type{4} = block_numbers(2:2:end,2:2:end);

for k=1:4
    A4{k} = A(type{k}(:));
end

%**** FUNCTION ****%
function [f_rev,P] = CoocReverse(f,T,order)
% Cooccurrences on reversed n-tuples, n=order
% f must be an array of (2*T+1)^order
% f_rev ... ccooccurrences for reversed order of input n-tuples, n=order
% P     ... permutation reversing the order of ccooccurrences
% Example: When f is the array of 4th order cooccurrences on [L,C,E,R]
%  Co = cooc([L(:),C(:),E(:),R(:)]',T);  f = Co(:);
%  then f_rev is the array of 4th order cooccurrences on [R,E,C,L]
%  The same array of cooccurrences is computed as 
% B = 2*T+1;  A = L(:)+T + B*(C(:)+T) + B^2*(E(:)+T) + B^3*(R(:)+T) + 1;  
% f = histc(A,1:B^4);
%  and also as
% fn = Cooc1split(D,4,'hor',T);   f = fn*prod(size(L));
% 
% EXAMPLE:
% B = 2*T+1;  A = L(:)+T + B*(C(:)+T) + B^2*(E(:)+T) + B^3*(R(:)+T) + 1;  f = histc(A,1:B^4);
% [f_rev,P] = CoocReverse(f,3,4);   
% B = 2*T+1;  A = R(:)+T + B*(E(:)+T) + B^2*(C(:)+T) + B^3*(L(:)+T) + 1;  ff = histc(A,1:B^4);
% all(f_rev==ff)
% ans =
%      1

B = (2*T+1);
if length(f)~= B^order, error('Inconsistent feature vector length with T and order'), end

f_rev = zeros(length(f),1);    P = zeros(length(f),1);

% for i=0:B^order-1
%     ntuple = dec2n(i,B,order) - T;          % shift the values to be in [-T,T] 
%     bin = BinNumber(ntuple+T,B,order);      % shift the values to be non-negative  
%     bin_rev = BinNumber(ntuple(end:-1:1)+T,B,order);
%     f_rev(bin_rev) = f(bin);
%     if nargout==2, 
%         P(bin_rev) = bin;
%     end
% end

for i=0:B^order-1
    ntuple = dec2n(i,B,order);            % values in [0,B-1] 
    bin = BinNumber(ntuple,B,order);  
    bin_rev = BinNumber(ntuple(end:-1:1),B,order);
    f_rev(bin_rev) = f(bin);
    if nargout==2, 
        P(bin_rev) = bin;
    end
end

%**** FUNCTION ****%
function bin = BinNumber(ntuple,B,order)

B_nary = B.^(0:order-1);        % B-nary numbers
% example:    order=4,  B_nary = [1,B,B^2,B^3];   
bin = dot(ntuple,B_nary)+1;     % values range from 0 to B-1


%**** FUNCTION ****%
function ntuple = dec2n(m,B,order)
% convert decimal number m to B-ary number
% this is the inverse function of BinNumber
% ntuple    ... B-ary number

ntuple = zeros(1,order);
for j=1:order
    ntuple(j) = mod(m,B); 
    m = floor(m/B);
end

%**** FUNCTION ****%
function Y = Quant(X,q,T)
% Quantization routine
% X ... variable to be quantized/truncated
% T ... threshold
% q ... quantization step (with type = 'scalar') or a vector of increasing
% non-negative integers outlining the quantization process.
% Y ... quantized/truncated variable
% Example 0: when q is a positive scalar, Y = trunc(round(X/q),T)
% Example 1: q = [0 1 2 3] quantizes 0 to 0, 1 to 1, 2 to 2, [3,Inf) to 3,
% (-Inf,-3] to -3, -2 to -2, -1 to -1. It is equivalent to Quant(.,3,1).
% Example 2: q = [0 2 4 5] quantizes 0 to 0, {1,2} to 1, {3,4} to 2,
% [5,Inf) to 3, and similarly with the negatives.
% Example 3: q = [1 2] quantizes {-1,0,1} to 0, [2,Inf) to 1, (-Inf,-2] to -1.
% Example 4: q = [1 3 7 15 16] quantizes {-1,0,1} to 0, {2,3} to 1, {4,5,6,7}
% to 2, {8,9,10,11,12,13,14,15} to 3, [16,Inf) to 4, and similarly the
% negatives.

if numel(q) == 1
    if q > 0, Y = trunc(round(X/q),T);
    else fprintf('*** ERROR: Attempt to quantize with non-positive step. ***\n'),end
else
    q = round(q);   % Making sure the vector q is made of integers
    if min(q(2:end)-q(1:end-1)) <= 0
        fprintf('*** ERROR: quantization vector not strictly increasing. ***\n')
    end
    if min(q) < 0, fprintf('*** ERROR: Attempt to quantize with negative step. ***\n'),end
    
    T = q(end);   % The last value determines the truncation threshold
    v = zeros(1,2*T+1);   % value-substitution vector
    Y = trunc(X,T)+T+1;   % Truncated X and shifted to positive values
    if q(1) == 0
        v(T+1) = 0; z = 1; ind = T+2;
        for i = 2 : numel(q)
            v(ind:ind+q(i)-q(i-1)-1) = z;
            ind = ind+q(i)-q(i-1);
            z = z+1;
        end
        v(1:T) = -v(end:-1:T+2);
    else
        v(T+1-q(1):T+1+q(1)) = 0; z = 1; ind = T+2+q(1);
        for i = 2 : numel(q)
            v(ind:ind+q(i)-q(i-1)-1) = z;
            ind = ind+q(i)-q(i-1);
            z = z+1;
        end
        v(1:T-q(1)) = -v(end:-1:T+2+q(1));
    end
    Y = v(Y);   % The actual quantization :)
end

%**** FUNCTION ****%
function Z = trunc(X,T)
% Truncation to [-T,T]
Z = X;
Z(Z > T)  =  T;
Z(Z < -T) = -T;

%**** FUNCTION ****%
function f_merged = Symm2split(f,T,order,type)
% length(f) is twice more due to splitting residuals to green and red/blue
%    pixels
le = length(f);
switch type
    case 'minmax'
        f1 = [f(1:le/4),     f(le/2+1:3/4*le)];   
        f2 = [f(le/4+1:le/2),f(3/4*le+1:end)];  
    case 'minmaxdir'
        f1 = [f(1:le/4),     f(le/2+1:3/4*le)];   
        f2 = [f(le/4+1:le/2),f(3/4*le+1:end)];  
    otherwise
        f1 = f(1:le/2);   f2 = f(le/2+1:end);  
end
f_merged = [Symm2(f1,T,order,type),Symm2(f2,T,order,type)];

%**** FUNCTION ****%
function [f_merged,EquivPair] = Symm2(f,T,order,type)
% Symm2.m  executes the symmetrizations in the same order as those in SRM,
%   SRMQ1, SRMQ1C, and SRMQ1CFA  (Outputs from Symm.m are in a different
%   order, i.e. permuted)  See the test Symm2T.m for an experimental proof.
% Merging cooccurrences of n-tuples, n=order based on symmetries
% type  ... 'dir', 'sign', 'dirsign', 'minmax', 'minmaxdir', 'signundo'
% f must be an array of length (2*T+1)^order or in case of 'minmax' and
%   'minmaxdir' two cooccurrence arrays are concatenated to twice this length 
% f_merged ... merged ccooccurrences 
%   Also see CoocReverse.m
%
% EXAMPLE:
% F_merged = Symm(F,2,4,'dirsign');   
%   is the same as 
% F_merged = symm1(F,2,4);
%
% F = Symm([Min Max],T,order,'minmaxdir');  
%     is the same as 
% F1 = Symm([Min Max],T,order,'minmax');  F = Symm(F1,T,order,'dir'); 

B = (2*T+1);   dim = B^order;
if nargin<4, error('The 4th input parameter ''type'' is required.'), end
if strcmp(type(end-2:end),'ndo'), dim = (B^order+1)/2; end 
if strcmp(type(1:3),'min'), dim = 2*B^order; end 
if length(f)~= dim, error('Inconsistent feature vector length with T and order'), end

f_merged = zeros(1,ceil(length(f)/2)); 
switch type
    case 'dirsign'
        EquivPair = zeros(length(f_merged),4);
    case 'minmaxdir'
        EquivPair = zeros(length(f_merged),4);
    otherwise
        EquivPair = zeros(length(f_merged),2);
end
m=0;  
switch type
case 'dir'
  % all zeros taken at first like at symm1.m    
  for i=0:B^order-1     % simple and correct
    ntuple = dec2n(i,B,order) - T;          % shift the values to be in [-T,T] 
    bin = BinNumber(ntuple+T,B,order);      % shift the values to be non-negative  
      % the only type specific two lines of code:
      bin_dire = BinNumber(ntuple(end:-1:1)+T,B,order); % reversed order in directional symmetry
      equivalency = [bin,bin_dire]; 
    if bin==min(equivalency)
        m=m+1;
        EquivPair(m,:) = equivalency;
        equivalency = unique(equivalency); 
        f_merged(m) = sum(f(equivalency));
    end
  end
case 'sign'
  done = zeros(1,B^order);
  for i=[(B^order-1)/2, 0:(B^order-1)/2-1,(B^order-1)/2+1:(B^order-1)]  
    ntuple = dec2nReversed(i,B,order) - T;  % shift the values to be in [-T,T] 
    bin = BinNumber(ntuple+T,B,order);      % shift the values to be non-negative  
      bin_sign = BinNumber(-ntuple+T,B,order);       
      equivalency = [bin,bin_sign]; 
    if ~done(bin)
        m=m+1;
        EquivPair(m,:) = equivalency;
        equivalency = unique(equivalency);  done(equivalency) = 1;
        f_merged(m) = sum(f(equivalency));
    end
  end
case 'dirsign'
  for i=[(B^order-1)/2, 0:(B^order-1)/2-1,(B^order-1)/2+1:(B^order-1)]  
    ntuple = dec2n(i,B,order) - T;          % shift the values to be in [-T,T] 
    bin = BinNumber(ntuple+T,B,order);      % shift the values to be non-negative  
      bin_dire = BinNumber(ntuple(end:-1:1)+T,B,order); % reversed order in directional symmetry
      bin_sign = BinNumber(-ntuple+T,B,order);       
      bin_dirsign  = BinNumber(-ntuple(end:-1:1)+T,B,order);
      equivalency = [bin,bin_sign,bin_dire,bin_dirsign];
    if bin==min(equivalency)
        m=m+1;
        EquivPair(m,:) = equivalency;
        equivalency = unique(equivalency);
        f_merged(m) = sum(f(equivalency));
    end
  end
case 'minmax'
  for i=0:B^order-1 
    ntuple = dec2n(i,B,order) - T;          % shift the values to be in [-T,T] 
    bin = BinNumber(ntuple+T,B,order);      % shift the values to be non-negative  
      bin_sign = BinNumber(-ntuple+T,B,order);       
    % two cooccurrence arrays are concatenated
      equivalency = [bin,bin_sign+B^order];
    if bin==min(equivalency)
        m=m+1;
        EquivPair(m,:) = equivalency;
        equivalency = unique(equivalency);
        f_merged(m) = sum(f(equivalency));
    end
  end
case 'minmaxdir'                
    fmm = Symm2(f,T,order,'minmax');
    [f_merged,EquivPair] = Symm2(fmm,T,order,'dir');  m = size(EquivPair,1);
case 'signundo'
  done = zeros(1,B^order);   k=1;
  for i=[(B^order-1)/2, 0:(B^order-1)/2-1,(B^order-1)/2+1:(B^order-1)]  
    ntuple = dec2nReversed(i,B,order) - T;  % shift the values to be in [-T,T] 
    bin = BinNumber(ntuple+T,B,order);      % shift the values to be non-negative  
      bin_sign = BinNumber(-ntuple+T,B,order);       
      equivalency = [bin,bin_sign]; 
    if ~done(bin)
         EquivPair(k,:) = equivalency;
         equivalency = unique(equivalency);  done(equivalency) = 1;
%         f_merged(m) = sum(f(equivalency));  m=m+1;
         f_merged(max(equivalency)) = f(k);  k=k+1;  % put the sum into one of the bins it came from
         m = max(m,max(equivalency));
    end
  end  
end
f_merged(m+1:end) = [];  EquivPair(m+1:end,:) = [];


function ntuple = dec2nReversed(m,B,order)
% convert decimal number m to B-ary number
% this is the inverse function of BinNumber
% ntuple    ... B-ary number
% version2: reversed order

ntuple = zeros(1,order);
for j=order:-1:1
    ntuple(j) = mod(m,B); 
    m = floor(m/B);
end


%**** FUNCTION ****%
function D = Residual(X,order,type)
% Computes the noise residual of a given type and order from MxN image X.
% residual order \in {1,2,3,4,5,6}
% type \in {hor,ver,diag,mdiag,KB,edge-h,edge-v,edge-d,edge-m}
% The resulting residual is an (M-b)x(N-b) array of the specified order,
% where b = ceil(order/2). This cropping is little more than it needs to 
% be to make sure all the residuals are easily "synchronized".
% !!!!!!!!!!!!! Use order = 2 with KB and all edge residuals !!!!!!!!!!!!!

[M N] = size(X);
I = 1+ceil(order/2) : M-ceil(order/2);
J = 1+ceil(order/2) : N-ceil(order/2);

switch type
    case 'hor'
        switch order
            case 1, D = - X(I,J) + X(I,J+1);
            case 2, D = X(I,J-1) - 2*X(I,J) + X(I,J+1);
            case 3, D = X(I,J-1) - 3*X(I,J) + 3*X(I,J+1) - X(I,J+2);
            case 4, D = -X(I,J-2) + 4*X(I,J-1) - 6*X(I,J) + 4*X(I,J+1) - X(I,J+2);
            case 5, D = -X(I,J-2) + 5*X(I,J-1) - 10*X(I,J) + 10*X(I,J+1) - 5*X(I,J+2) + X(I,J+3);
            case 6, D = X(I,J-3) - 6*X(I,J-2) + 15*X(I,J-1) - 20*X(I,J) + 15*X(I,J+1) - 6*X(I,J+2) + X(I,J+3);
        end
    case 'ver'
        switch order
            case 1, D = - X(I,J) + X(I+1,J);
            case 2, D = X(I-1,J) - 2*X(I,J) + X(I+1,J);
            case 3, D = X(I-1,J) - 3*X(I,J) + 3*X(I+1,J) - X(I+2,J);
            case 4, D = -X(I-2,J) + 4*X(I-1,J) - 6*X(I,J) + 4*X(I+1,J) - X(I+2,J);
            case 5, D = -X(I-2,J) + 5*X(I-1,J) - 10*X(I,J) + 10*X(I+1,J) - 5*X(I+2,J) + X(I+3,J);
            case 6, D = X(I-3,J) - 6*X(I-2,J) + 15*X(I-1,J) - 20*X(I,J) + 15*X(I+1,J) - 6*X(I+2,J) + X(I+3,J);
        end
    case 'diag'
        switch order
            case 1, D = - X(I,J) + X(I+1,J+1);
            case 2, D = X(I-1,J-1) - 2*X(I,J) + X(I+1,J+1);
            case 3, D = X(I-1,J-1) - 3*X(I,J) + 3*X(I+1,J+1) - X(I+2,J+2);
            case 4, D = -X(I-2,J-2) + 4*X(I-1,J-1) - 6*X(I,J) + 4*X(I+1,J+1) - X(I+2,J+2);
            case 5, D = -X(I-2,J-2) + 5*X(I-1,J-1) - 10*X(I,J) + 10*X(I+1,J+1) - 5*X(I+2,J+2) + X(I+3,J+3);
            case 6, D = X(I-3,J-3) - 6*X(I-2,J-2) + 15*X(I-1,J-1) - 20*X(I,J) + 15*X(I+1,J+1) - 6*X(I+2,J+2) + X(I+3,J+3);
        end
    case 'mdiag'
        switch order
            case 1, D = - X(I,J) + X(I-1,J+1);
            case 2, D = X(I-1,J+1) - 2*X(I,J) + X(I+1,J-1);
            case 3, D = X(I-1,J+1) - 3*X(I,J) + 3*X(I+1,J-1) - X(I+2,J-2);
            case 4, D = -X(I-2,J+2) + 4*X(I-1,J+1) - 6*X(I,J) + 4*X(I+1,J-1) - X(I+2,J-2);
            case 5, D = -X(I-2,J+2) + 5*X(I-1,J+1) - 10*X(I,J) + 10*X(I+1,J-1) - 5*X(I+2,J-2) + X(I+3,J-3);
            case 6, D = X(I-3,J+3) - 6*X(I-2,J+2) + 15*X(I-1,J+1) - 20*X(I,J) + 15*X(I+1,J-1) - 6*X(I+2,J-2) + X(I+3,J-3);
        end
    case 'KB'
        D = -X(I-1,J-1) + 2*X(I-1,J) - X(I-1,J+1) + 2*X(I,J-1) - 4*X(I,J) + 2*X(I,J+1) - X(I+1,J-1) + 2*X(I+1,J) - X(I+1,J+1);
    case 'edge-h'
        Du = 2*X(I-1,J) + 2*X(I,J-1) + 2*X(I,J+1) - X(I-1,J-1) - X(I-1,J+1) - 4*X(I,J);   %   -1  2 -1
        Db = 2*X(I+1,J) + 2*X(I,J-1) + 2*X(I,J+1) - X(I+1,J-1) - X(I+1,J+1) - 4*X(I,J);   %    2  C  2    +  flipped vertically
        D = [Du,Db];
    case 'edge-v'
        Dl = 2*X(I,J-1) + 2*X(I-1,J) + 2*X(I+1,J) - X(I-1,J-1) - X(I+1,J-1) - 4*X(I,J);   %   -1  2
        Dr = 2*X(I,J+1) + 2*X(I-1,J) + 2*X(I+1,J) - X(I-1,J+1) - X(I+1,J+1) - 4*X(I,J);   %    2  C       +  flipped horizontally
        D = [Dl,Dr];                                                                      %   -1  2
    case 'edge-m'
        Dlu = 2*X(I,J-1) + 2*X(I-1,J) - X(I-1,J-1) - X(I+1,J-1) - X(I-1,J+1) - X(I,J); %      -1  2 -1
        Drb = 2*X(I,J+1) + 2*X(I+1,J) - X(I+1,J+1) - X(I+1,J-1) - X(I-1,J+1) - X(I,J); %       2  C       +  flipped mdiag
        D = [Dlu,Drb];                                                                 %      -1
    case 'edge-d'
        Dru = 2*X(I-1,J) + 2*X(I,J+1) - X(I-1,J+1) - X(I-1,J-1) - X(I+1,J+1) - X(I,J); %      -1  2 -1
        Dlb = 2*X(I,J-1) + 2*X(I+1,J) - X(I+1,J-1) - X(I+1,J+1) - X(I-1,J-1) - X(I,J); %          C  2    +  flipped diag
        D = [Dru,Dlb];                                                                 %            -1
    case 'KV'
        D = 8*X(I-1,J) + 8*X(I+1,J) + 8*X(I,J-1) + 8*X(I,J+1);
        D = D - 6*X(I-1,J+1) - 6*X(I-1,J-1) - 6*X(I+1,J-1) - 6*X(I+1,J+1);
        D = D - 2*X(I-2,J) - 2*X(I+2,J) - 2*X(I,J+2) - 2*X(I,J-2);
        D = D + 2*X(I-1,J-2) + 2*X(I-2,J-1) + 2*X(I-2,J+1) + 2*X(I-1,J+2) + 2*X(I+1,J+2) + 2*X(I+2,J+1) + 2*X(I+2,J-1) + 2*X(I+1,J-2);
        D = D - X(I-2,J-2) - X(I-2,J+2) - X(I+2,J-2) - X(I+2,J+2) - 12*X(I,J);
end
