%clear console; Remove all variables from memory; Close all figures
clc;clearvars;close all;

% MIMO system

m = 4;
n = 4;
l = 16;

alphabet = [ exp(1i*pi/4); exp(1i*3/4*pi); exp(1i*5/4*pi); exp(1i*7/4*pi) ];

Number_Of_CR = 2000;%1000
Number_Of_increasing_PTX = 9;%7

SumMSEinit = zeros(1,Number_Of_increasing_PTX);
SumMSE_LD = zeros(1,Number_Of_increasing_PTX);
SumMSE_SD = zeros(1,Number_Of_increasing_PTX);
SumMSE_VB = zeros(1,Number_Of_increasing_PTX);

Sumiter_ld = zeros(1,Number_Of_increasing_PTX);
Sumiter_sd = zeros(1,Number_Of_increasing_PTX);
Sumiter_vb = zeros(1,Number_Of_increasing_PTX);

for CR = 1:Number_Of_CR

    CR
    h0 = crandn(m, n);
    h1 = crandn(n, l);
    h2 = crandn(m, l);
    vinit = alphabet(round(4*rand(l, 1)+0.5)); %understood

    ptxdB = -30;%-30
    ptx = db2mag(ptxdB);%must be very low, have to reformulate ptx.
    for i =  1:Number_Of_increasing_PTX
        ptxdB;

        ptx = ptx/n;

        htotinit = (h0+h2*diag(vinit)*h1.')*sqrt(ptx); %understood %note: transpose(h1) == h1.'
        ginit = (eye(n)+htotinit'*htotinit)^-1*htotinit';
        mseinit = real(trace((eye(n)-ginit*htotinit)'*(eye(n)-ginit*htotinit)+ginit*ginit'));%mseinit must equal checkmseinit
        %checkmseinit = trace((eye(n)+htotinit'*htotinit)^-1)
        %checkmseinit2 = (vec(eye(n)-ginit*htotinit))'*vec(eye(n)-ginit*htotinit)+trace(ginit*ginit')

        %%%%%%%%%%Check Eqns%%%%%%%%%%
        b_sdinit = vec(eye(n)-sqrt(ptx)*ginit*h0);
        a_sdinit = sqrt(ptx)*kr(h1,ginit*h2);
        Rewritten_mseinit = (b_sdinit-a_sdinit*vinit)'*(b_sdinit-a_sdinit*vinit)+trace(ginit*ginit');
        if norm(Rewritten_mseinit-mseinit) < 1e-8
            %disp("mse eqn for MLD is written correctly")
        else
            error("mse eqn for MLD is not written correctly")
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%% no IRS %%%%%%%%%%%%%%%%%%%%%%%%%

        SumMSEinit(i) = SumMSEinit(i) + mseinit;

        %%%%%%%%%%%%%%%%%%%%%% V-BLAST Detection: Pvb*ECMvb_init*transpose(Pvb)=Lvb*Dvb*Lvb' %%%%%%%%%%%%%%%%%%%%%%%%%
        gvbold = ginit;
        %vvbold = vinit;
        msevbold = mseinit;
        iter_vb = 1;
        while(1)
            b_vb = vec(eye(n)-sqrt(ptx)*gvbold*h0);
            a_vb = sqrt(ptx)*kr(h1,gvbold*h2);
            %B = a_vb'*a_vb+(trace(gvbold*gvbold')/l)*eye(l);
            %vblue = B^-1*a_vb'*b_vb;
            %initialize Permutation & Diagonal matrices
            Pvb = eye(l);
            Dvb = zeros(l);
            B = a_vb'*a_vb+(trace(gvbold*gvbold')/l)*eye(l);
            ECMvb = (B)^-1;        ECMvb_init=(B)^-1;
            Identity = eye(l);
            for DS = 1:l %DS = DataStream %If there is 4 Transmit Antennas, there is 4 DataStreams getting transmitted.
                DS;
                ECMvb; %check
                MSE = diag(ECMvb);
                MSE = MSE(DS:l,:);
                minMSE = min(MSE);
                index = find(MSE==minMSE);
                index = index + (DS-1);%correct
                Pi = Identity;
                Pi([DS index],:) = Pi([index DS],:);%correct %DS-th and index rows are exchanged
                Pvb = Pi*Pvb;
                ECMvb = Pi * ECMvb * transpose(Pi);
                Dvb(DS,DS) = ECMvb(DS,DS);
                ECMvb(DS:l,DS) = ECMvb(DS:l,DS)/Dvb(DS,DS);
                for j = DS+1:l
                    ECMvb(j:l,j) = ECMvb(j:l,j) - ECMvb(j:l,DS) * conj(ECMvb(j,DS)) * Dvb(DS,DS);
                    ECMvb(j,j:l) = ECMvb(j:l,j)';
                end
            end
            Dvb;
            Lvb=tril(ECMvb);
            LHS=Pvb*ECMvb_init*transpose(Pvb);
            RHS=Lvb*Dvb*Lvb';
            if norm(LHS-RHS)<1e-8
                %disp("CF is done correctly")
            else
                error("CF is not done correctly")
            end
            %Forward & Backward Filters
            Bfiltervb = (Lvb^-1)';
            Ffiltervb = ( Dvb * Lvb' * Pvb * B' )';
            %Use Filters to find Svblast
            Vvb_p = zeros(l,1);
            Vvb_p_q = zeros(l,1);
            bfilter_vb = Bfiltervb' - eye(l);
            for roww = 1:l
                Vvb_p(roww,1) = Ffiltervb(:,roww)' * (a_vb')*b_vb - bfilter_vb(roww,:) * Vvb_p_q;
                Vvb_p_q(roww,1) = Quantizer(Vvb_p(roww,1));
            end
            Vvb_p;
            Vvb_p_q;
            Vvb_q = transpose(Pvb) * Vvb_p_q;
            Vvb = Vvb_q; %Checked: correct
            Vvb;
            Vvbnew = Vvb;
            htotvbnew = (h0+h2*diag(Vvbnew)*h1.')*sqrt(ptx);
            gvbnew = (eye(n)+htotvbnew'*htotvbnew)^-1*htotvbnew';
            msevbnew = trace((eye(n)+htotvbnew'*htotvbnew)^-1);
            if msevbnew >= msevbold
                %disp('MSE does not decrease')
                break;%understood
            end
            %vold = vldnew;
            gvbold = gvbnew;
            msevbold = msevbnew;
            iter_vb = iter_vb+1;
        end
        mse_vb = msevbold;
        SumMSE_VB(i) = SumMSE_VB(i) + mse_vb;



        %%%%%%%%%%%%%%%%%%%%%% Sphere Decoder %%%%%%%%%%%%%%%%%%%%%%%%%
        gsdold = ginit;
        %vsdold = vinit;
        msesdold = mseinit;
        iter_sd = 1;
        while (1)
            b_sd = vec(eye(n)-sqrt(ptx)*gsdold*h0);
            a_sd = sqrt(ptx)*kr(h1,gsdold*h2);
            global Lsd;
            global Dsd;
            global est_X;
            global V;
            global r;
            global Vsd;
            %initialize Permutation & Diagonal matrices
            Dsd = zeros(l);
            B = a_sd'*a_sd+(trace(gsdold*gsdold')/l)*eye(l);       Binit = a_sd'*a_sd+(trace(gsdold*gsdold')/l)*eye(l);%for checking purpose: norm(P*ahainit*transpose(P)-L'*D*L)
            Identity = eye(l);
            for vs = l:-1:1 %DS = DataStream %If there is 4 Transmit Antennas, there is 4 DataStreams getting transmitted.
                vs;
                B; %check
                Dsd(vs,vs) = B(vs,vs);
                B(1:vs,vs) = B(1:vs,vs)/Dsd(vs,vs);
                for j = 1:(vs-1)
                    B(1:j,j) = B(1:j,j) - B(1:j,vs) * conj(B(j,vs)) * Dsd(vs,vs);
                    B(j,1:j) = B(1:j,j)';
                end
            end
            Dsd;
            Lsd = (triu(B))';
            CFeqn = Lsd'*Dsd*Lsd;
            if norm(Binit-CFeqn)<1e-8 %eqninit(1,:) CFeqn(1,:)
                %disp("CF is done correctly")
            else
                error("CF is not done correctly")
            end
            est_X = Dsd^-1*(Lsd')^-1*a_sd'*b_sd; %X is consistently updated

            %Possible_Signal at each level
            V = [exp(1i*(pi/4)) exp(1i*(3.*pi/4)) exp(1i*(5.*pi/4)) exp(1i*(7.*pi/4))];
            Vsd = zeros(l,1);
            for ii = 1:1 %No. of symbols = 1
                r = inf;
                iii = 1; %start with node = 1
                sumMSE = 0;
                SiPath = zeros(l,1);
                SDnode(ii, iii, sumMSE, SiPath,l);
            end
            Vsd;
            vsdnew = Vsd;
            htotsdnew = (h0+h2*diag(vsdnew)*h1.')*sqrt(ptx);
            gsdnew = (eye(n)+htotsdnew'*htotsdnew)^-1*htotsdnew';
            %msesdnew = real(trace((eye(n)-gsdnew*htotsdnew)'*(eye(n)-gsdnew*htotsdnew)+gsdnew*gsdnew'));
            msesdnew = trace((eye(n)+htotsdnew'*htotsdnew)^-1);
            
            if msesdnew >= msesdold
                %disp('MSE does not decrease')
                break;%understood
            end
            %vold = vsdnew;%understood
            gsdold = gsdnew;%understood
            msesdold = msesdnew;%understood
            iter_sd = iter_sd+1;
        end
        mse_sd = msesdold;
        SumMSE_SD(i) = SumMSE_SD(i) + mse_sd;


        %%%%%%%%%%%%%%%%%%%% Linear Detection %%%%%%%%%%%%%%%%%%%%%%%
        gldold = ginit;
        %vldold = vinit;
        mseldold = mseinit;
        iter_ld = 1;
        while (1)
            b_ld = vec(eye(n)-sqrt(ptx)*gldold*h0);
            a_ld = sqrt(ptx)*kr(h1,gldold*h2);
            %initialize Permutation & Diagonal matrices
            Dld = zeros(l);
            B = a_ld'*a_ld+(trace(gldold*gldold')/l)*eye(l); 
            vblue = B^-1*a_ld'*b_ld;
            vlin = sqrt(2)*((real(vblue) > 0)*1-0.5+1i*((imag(vblue) > 0)*1-0.5)); %Quantization check: [vblue vldnew].'
            vldnew = vlin;
            htotldnew = (h0+h2*diag(vldnew)*h1.')*sqrt(ptx);
            gldnew = (eye(n)+htotldnew'*htotldnew)^-1*htotldnew';
            %mseldnew = real(trace((eye(n)-gldnew*htotldnew)'*(eye(n)-gldnew*htotldnew)+gldnew*gldnew'));
            mseldnew = trace((eye(n)+htotldnew'*htotldnew)^-1);
            if mseldnew >= mseldold
                %disp('MSE does not decrease')
                break;%understood
            end
            %vold = vldnew;
            gldold = gldnew;
            mseldold = mseldnew;
            iter_ld = iter_ld+1;
        end
        mse_ld = mseldold;
        SumMSE_LD(i) = SumMSE_LD(i) + mse_ld;

        %         iter_ld;
        %         iter_sd;
        %
        %         %Note: if Sumiter = CR, v didn't change at all for that ptxdB value
        %         Sumiter_ld(i) = Sumiter_ld(i) + iter_ld;
        %         Sumiter_sd(i) = Sumiter_sd(i) + iter_sd;



        ptxdB = ptxdB + 10;
        ptx = db2mag(ptxdB);
    end
end


%PLOTS
ptxdB_plot = zeros(1, Number_Of_increasing_PTX); %Coordinates
AVE_MSE_LD = zeros(1, Number_Of_increasing_PTX);
AVE_MSE_SD = zeros(1, Number_Of_increasing_PTX);
AVE_MSE_init = zeros(1, Number_Of_increasing_PTX);
AVE_MSE_VB = zeros(1, Number_Of_increasing_PTX);

ptxdB = -30;%-30
ptx = db2mag(ptxdB);
for i =  1:Number_Of_increasing_PTX

    ptxdB;
    ptx;

    ptxdB_plot(i) = ptxdB;
    AVE_MSE_LD(i) = SumMSE_LD(i) / Number_Of_CR;
    AVE_MSE_SD(i) = SumMSE_SD(i) / Number_Of_CR;
    AVE_MSE_init(i) = SumMSEinit(i) / Number_Of_CR;
    AVE_MSE_VB(i) = SumMSE_VB(i) / Number_Of_CR;

    ptxdB = ptxdB + 10;
    ptx = db2mag(ptxdB);

end

% AVE_MSE_init
% AVE_MSE_LD
% AVE_MSE_SD


figure(1)
semilogy(ptxdB_plot, real(AVE_MSE_init),'color','blue') % line and marker
hold on
semilogy(ptxdB_plot, real(AVE_MSE_LD),'color','green') % line and marker
hold on
semilogy(ptxdB_plot, real(AVE_MSE_VB),'color','black') % line and marker
hold on
semilogy(ptxdB_plot, real(AVE_MSE_SD),'color','red') % line and marker
hold on

legend("no irs","LD","VB","SD")
% legend("no irs","VB")

title("MSE(PTX)");
xlabel('PTX(dB-scale)'); %want it to be dB
ylabel('MSE(log-scale)')


function SDnode(ii, iii, sumMSE, SiPath, TxRx)    %ii is impt to pass the value of the current symbol No. pass summation of Mu to current node
global r;
global Vsd;
global Dsd;
global V;
global est_X;
Viewed = zeros(4,1);
CheckMinMSEIndex = zeros(4,1);
sumLsF = sumLs(iii, SiPath);
while (min(Viewed) == 0)     %not all 4 neighbours path is check yet
    %iii
    for k=1:4
        if Viewed(k) == 0
            CheckMinMSEIndex(k) = real( Dsd(iii,iii) * norm( est_X(iii,ii) - V(k) - sumLsF)^2 );%Find the 4 neighbours MSE value
            %MSE3(1,n) = real( Dsd(3,3) * norm( est_X(3,ii) - S3(n) - Lsd(3,1)*minMSE_S1 - Lsd(3,2)*minMSE_S2 )^2 );
        else
            CheckMinMSEIndex(k) = inf;
        end
    end
    MinMSEIndex = find(CheckMinMSEIndex == min(CheckMinMSEIndex));%Find the index of the current node MSE\
    %logic = (norm(est_X(iii,ii) - V(MinMSEIndex) - sumLsF))^2 <= (r^2 - sumMSE)/D(iii,iii) %suppose to be logical 1 until iii=l
    %LHS = (norm(est_X(iii,ii) - V(MinMSEIndex) - sumLsF))^2
    %RHS = (r^2 - sumMSE)/D(iii,iii)
    %sizeofD = size(D)
    %sizeofest_X = size(est_X)
    %sizeofV = size(V)
    %sizeofr = size(r)
    %sizeofVsd = size(Vsd)
    %sizeofsumLsF = size(sumLsF)
    if (norm(est_X(iii,ii) - V(MinMSEIndex) - sumLsF))^2 <= (r^2 - sumMSE)/Dsd(iii,iii)
        NextSiPath = SiPath;
        NextSiPath(iii) = V(MinMSEIndex);
        if iii ~= TxRx
            Viewed(MinMSEIndex) = 1;
            SDnode(ii, iii+1, sumMSE + min(CheckMinMSEIndex), NextSiPath, TxRx);
        else
            r = sqrt(sumMSE + min(CheckMinMSEIndex));
            Vsd(:,ii) = NextSiPath;
            break;
        end
        %      elseif (norm(est_X(iii,ii) - S(MinMSEIndex) - sumLsF))^2 > (r^2 - sumMSE)/Dsd(iii,iii)
    else
        break;
    end
end
end

function sumLsF = sumLs(iii, SiPath)
global Lsd;
tempSum = 0;
for j = 1:iii-1
    tempSum = tempSum + Lsd(iii,j) * SiPath(j); %SiPath is the minMSE_Sj
end
sumLsF = tempSum;
end


function R = vec(A) %Vectorization(mathematics)

R = reshape(A,numel(A), 1);

end

function AB = kr(A,B)%khatri-rao product [Column-wise Kronecker product]
[I,F]=size(A);
[J,F1]=size(B);
if F~=F1
    error(' Error in kr.m - The matrices must have the same number of columns')
end
AB=zeros(I*J,F);
for f=1:F
    ab=B(:,f)*A(:,f).';
    AB(:,f)=ab(:);
end
end

function QuantizerF = Quantizer(element)
elementRe = real(element);
elementIm = imag(element);
if elementRe>0 && elementIm>0
    element_q = exp(1i*(pi/4));
elseif elementRe<0 && elementIm>0
    element_q = exp(1i*(3*pi/4));
elseif elementRe<0 && elementIm<0
    element_q = exp(1i*(5*pi/4));
else
    element_q = exp(1i*(7*pi/4));
end
QuantizerF = element_q;
end
