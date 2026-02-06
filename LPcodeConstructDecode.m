%% Creating and decoding a lifted-product QLDPC code. 
% Creates and simulates the min-sum decoding of a LP-QLDPC code
% Starting from doing the same for a QC-LDPC code. 
%% Starting with the classial (3,5)-quasi-cyclic code with lift size 31:
l=31;
x=circshift(eye(l),1);
protoBpowers=[1 2 4 8 16;5 10 20 9 18;25 19 7 14 28];
%% Construction of the parity check matrix from the base protograph via lifting
global HB HBcopy HBheight HBwidth rowspan colspan genmat cn_neighbors vn_neighbors;
HB=ArrayWithHandle();
HB.data=sparse(cell2mat(arrayfun(@(int) x^int,protoBpowers,'UniformOutput',false)));
HBcopy=HB.data;
HBheight=height(HB.data);HBwidth=width(HB.data);
rowspan=1:HBheight;colspan=1:HBwidth;
cn_neighbors=cellfun(@(row) colspan(row>0), mat2cell(HB.data',HBwidth,ones(1,HBheight)),'UniformOutput',false);
vn_neighbors=cellfun(@(col) rowspan(col>0), mat2cell(HB.data,HBheight,ones(1,HBwidth)),'UniformOutput',false);
genmat=reducePCmat();

%Plotting the generator matrix
% imagesc(HBcopy);
imagesc([1,155],[1,64],genmat.data);
% axis equal;
% colorbar;
paritychecks=mod(genmat.data*(HBcopy'),2);
fprintf("\nAny parity check violations exist: ");
disp(logical(full(any(paritychecks,"all"))));
%% Simulation of the decoding of the classical quasi-cyclic code
%syndrome_MSA_seq_vars_5(HB.data,vn_neighbors,cn_neighbors,randi([0,1],1,HBheight),log2(9),20);<--- Test of minsum function
pvals=logspace(log10(.05),log10(.2),10);
errorprobs=arrayfun(@(p) simulateMinSum(p,0,40),pvals);
% errorprob=simulateMinSum(.1,0,20);
disp(errorprobs);

%% Generating the LP-QLDPC code by lifting the CSS-PC matrix associated with B
%Unlifted stabilizer matrix
mb=height(protoBpowers);nb=width(protoBpowers);
BX=[kron(protoBpowers,eye(nb)), kron(eye(mb),(l-protoBpowers)')];
BZ=[kron(eye(nb),protoBpowers),kron((l-protoBpowers)',eye(mb))];
% Following line generates the CSS-partity check (stabilizer) matrix.
HB.data=sparse(cell2mat(arrayfun(@(int) (int==0)*zeros(l)+(int~=0)*x^int,[BX,zeros(nb*mb,nb^2+mb^2);zeros(nb*mb,nb^2+mb^2),BZ],'UniformOutput',false)));

HBheight=height(HB.data);HBwidth=width(HB.data);
rowspan=1:HBheight;colspan=1:HBwidth;
global omega;
omega=spdiags(1,[-HBwidth/2-1,HBwidth/2+1],HBwidth,HBwidth);
% omega=[spalloc(HBwidth/2,HBwidth/2,0),spdiags(1,0,HBwidth/2,HBwidth/2);spdiags(1,0,HBwidth/2,HBwidth/2),spalloc(HBwidth/2,HBwidth/2,0)];
%Uncomment and run below two lines to row-reduce the stabilizer matrix of the LPQLDPC code and obtain the # of indep stabilizers.
% quantumgenmat=reducePCmat(); % Note 'quantumgenmat' does not correspond to the stabilized code subspace of the Hilbert Space
% fprintf('Number of independent stabilizers: %d\n',sum(diag(HB.data))); %Number of logical qubits is then HBwidth/2 less the trace (140 below).

%% Simulation of Min-Sum decoding of the LP-QLDPC code (Need to run above section before every run of the sections under this title)
% First, generate the logical Pauli Operators, assuming Z-type Z-logicals. (they live in rowsp(H(BX))^\perp)
%% 
% 1) Take 140 rows of G(H(BX)). They are automatically symplectically-mutually orthogonal.
HB.data=HB.data(1:HBheight/2,1:HBwidth/2); %assigning HX to HB.data
HBheight=height(HB.data);HBwidth=width(HB.data);
rowspan=1:HBheight;colspan=1:HBwidth;
LZgenmat=reducePCmat();
LZs=[zeros(140,HBwidth),LZgenmat.data(1:140,:)]; %A generator matrix for logical Zs (in symplectic representation)
% Re-assign HB to the CSS matrix by running the previous section before running step 2.

%% 
% 2) Solve for the j^th X-logical vector by adding LZs as rows to H and augmenting H by a column of zeros that is 1 on the j^th added row.
% CSSHcopy=HB.data;
num_lbits=140;
num_Xconditions=height(LZs)+HBheight;

HB.data=[HB.data;LZs];
HBheight=height(HB.data);HBwidth=width(HB.data);
rowspan=1:HBheight;colspan=1:HBwidth;

HB.data=HB.data*omega;
[Xgens,rowops,colops]=reducePCmat();
rowrankHBLZ=sum(diag(HB.data));
reduced_augvecs=rowops*[zeros(HBheight,num_lbits);eye(num_lbits)];
LXs=(augvecs(1:rowrankHBLZ,:)+(Xgens(1,:)/colops).').'*colops;

% 3) Symplectically orthogonalize LXs using the LQ decomposition for binary symplectic matrices
L=zeros(num_lbits);L(1,1)=1;
for j=2:num_lbits
    L(j,:)=LXs(1:j-1,:)*omega*LXs(j,:).';
    LXs(j,:)=LXs(j,:)-L*LXs;
end
LZs=L*LZs; %Since L is symplectic, joint transformation of LXs and LZs preserve their mutual commutation relations.


%% 
%Alternatively, get X-logicals from Narayanan's file:
global LX;
LX=read_L('LP_Tanner155_lx.alist');
fprintf('Logicals are inconsistent with Z-type stabilizers: %d\n',full(any(mod(LX*HB.data(446:end,1055:end).',2),'all')));

%% 
% Error performance Simulation
HB.data=HB.data(1:475,1:1054);
HBheight=height(HB.data);HBwidth=width(HB.data);
rowspan=1:HBheight;colspan=1:HBwidth;
cn_neighbors=cellfun(@(row) colspan(row>0), mat2cell(HB.data.',HBwidth,ones(1,HBheight)),'UniformOutput',false);
vn_neighbors=cellfun(@(col) rowspan(col>0), mat2cell(HB.data,HBheight,ones(1,HBwidth)),'UniformOutput',false);
pvals=logspace(log10(.04),log10(.09),6);
errorprobs=arrayfun(@(p) SimulateMinSum_Zchannel(p,45),pvals);
%Plotting
plot(pvals,errorprobs,'-o','MarkerFaceColor','blue')
xlabel('Z-flip probability')
ylabel('Decoding error probability')
title('Error performance of the [[1054,140,20]]-LP-QLDPC code against Z-errors')
grid on

%% Helper functions
function errorprob = simulateMinSum(p,b,nsamples) %Simulate decoding of LDPC code 
% defined by HB, for transmission through a BSC with t.p. p and bias b.
    global HBcopy genmat vn_neighbors cn_neighbors;
    failcount=0;
    for i=1:nsamples
        codedmsgsamp=mod(randi([0,1],1,height(genmat.data))*genmat.data,2);
        channeloutsamp=arrayfun(@(x) BSCwithbias(x,p,b),codedmsgsamp);
        error_vec=mod(channeloutsamp+codedmsgsamp,2);
        syndrome=mod(channeloutsamp*(HBcopy'),2);
        % fprintf('No of channel errors: %d\n',sum(channeloutsamp~=codedmsgsamp));
        % disp(strcat('Syndrome vector: ',num2str(syndrome)));
        minsum_error_guess=syndrome_MSA_seq_vars_5(HBcopy,vn_neighbors,cn_neighbors,syndrome,log2((1-p)/p),40);
        % disp(strcat('coded message: ',num2str(codedmsgsamp)));
        % disp(strcat('errors vector: ',num2str(error_vec)));
        % disp(strcat('min sum guess: ',num2str(minsum_error_guess)));
        % disp(strcat('mismachvector: ',num2str(codedmsgsamp~=minsum_error_guess)));fprintf('\n');
        if any(error_vec~=minsum_error_guess)
            % disp('decode fail');
            failcount=failcount+1;
        end
    end
    errorprob=failcount/nsamples;
end

function errorprob = SimulateMinSum_Zchannel(p,nsamples)
    global HB LX vn_neighbors cn_neighbors; %omega
    failcount=0;
    for i=1:nsamples
        Z_error_samp=arrayfun(@(x) BSCwithbias(x,p,0),zeros(1,width(HB.data)));%[zeros(1,width(HB.data)/2),arrayfun(@(x) BSCwithbias(x,p,0),zeros(1,width(HB.data)/2))];
        % disp(strcat('Z_error_dimensions: ',num2str(size(Z_error_samp))));
        syndrome=mod(Z_error_samp*HB.data.',2);%(omega*(HB.data.')),2);
        % disp(strcat('Syndrome vector: ',num2str(syndrome)));
        minsum_error_guess=syndrome_MSA_seq_vars_5(HB.data,vn_neighbors,cn_neighbors,syndrome,log2((1-p)/p),30);
        % disp(strcat('errors vector: ',num2str(Z_eror_samp)));
        % disp(strcat('min sum guess: ',num2str(size(minsum_error_guess))));
        if any(mod(mod(minsum_error_guess+Z_error_samp,2)*LX.',2),'all')%*omega*[LX,zeros(140,width(HB.data)/2)].',2),'all')
            % disp('decode fail');
            failcount=failcount+1;
        end
    end
    errorprob=failcount/nsamples;
    fprintf('Evaluated error: %d\n',errorprob);
end

function y = BSCwithbias(x,p,b)
%Function simulating the i.i.d. BSC with transition probability p and bias
%b. Given a channel input, randomly generates the channel output y as 0 ->
%1 with probability p+b and 1->0 with probability p-b. Make sure p and b
%are between 0 and 1.
    if rand(1)>(1-p)+(-1).^x*b
        y=mod(x+1,2);
    else 
        y=x;
    end
end

function [genmat,rowOPSmatrix,colOPSmatrix] = reducePCmat() 
%returns the generator matrix of an LDPC code after
%copyless row-reduction (mod w) of the global parity check matrix HB
    global HB HBheight HBwidth rowspan colspan;
    rowOPSmatrix=spalloc(HBheight,HBheight,3*HBheight); %Preallocation for performance optimization
    rowOPSdiag=sub2ind(size(rowOPSmatrix),1:HBheight,1:HBheight);
    rowOPSmatrix(rowOPSdiag)=ones(1,HBheight); %right-multiplication by this applies the row operations
    colOPSmatrix=sparse(eye(HBwidth)); %left-multiplication by this reverses the column operations
    % disp("I'm being called");
    shortside=min(HBheight,HBwidth);
    colswaps=[];
    for r=1:shortside% lower triangulation
        lowercolsupport=rowspan([zeros(r-1,1);HB.data(r:end,r)]'>0);
        % disp(lowercolsupport);
        if ~isempty(lowercolsupport)
            if lowercolsupport(1)~=r %lower row-swap in case zero on rth diagonal entry
                fprintf('Swapping rows %d and %d\n\n',r,lowercolsupport(1));
                HB.data(r,:)=mod(HB.data(r,:)+HB.data(lowercolsupport(1),:),2);
                HB.data(lowercolsupport(1),:)=mod(HB.data(lowercolsupport(1),:)+HB.data(r,:),2);
                HB.data(r,:)=mod(HB.data(r,:)+HB.data(lowercolsupport(1),:),2);
                %repeat for rowOPSmatrix
                rowOPSmatrix(r,:)=mod(rowOPSmatrix(r,:)+rowOPSmatrix(lowercolsupport(1),:),2);
                rowOPSmatrix(lowercolsupport(1),:)=mod(rowOPSmatrix(lowercolsupport(1),:)+rowOPSmatrix(r,:),2);
                rowOPSmatrix(r,:)=mod(rowOPSmatrix(r,:)+rowOPSmatrix(lowercolsupport(1),:),2);
            end
        else
            % disp("Triggered empty lower column condition");
            % fprintf("Upper row at column %d: ",r);fprintf(strcat(num2str(HB.data(r,r+1:end)),'\n'));
            upperrowsupport=colspan([zeros(1,r),HB.data(r,r+1:end)]>0);
            % fprintf("Upper row support at column %d: ",r);fprintf(strcat(num2str(upperrowsupport),'\n'));
            fprintf(strcat("Upper row support is empty: ",num2str(isempty(upperrowsupport)),"\n"));
            if isempty(upperrowsupport) % column-swap in case zero on lower column
                %Send empty row to bottom of matrix 
                j=0;
                while ~any(HB.data(end-j,:)) && HBheight-j>r%Need to make sure bottom row is not all zeros.
                    j=j+1;
                end                    
                if HBheight-j>r
                    fprintf(strcat("Sending row ",num2str(r)," to bottom (swapping with row %d)\n"),HBheight-j);
                    HB.data(r,:)=HB.data(end-j,:);%mod(HB.data(r,:)+HB.data(end,:),2);
                    HB.data(end-j,:)=zeros(1,HBwidth);%mod(HB.data(end,:)+HB.data(r,:),2);
                    % HB.data(r,:)=mod(HB.data(r,:)+HB.data(end,:),2);
                    %repeat for rowOPSmatrix
                    rowOPSmatrix(r,:)=mod(rowOPSmatrix(r,:)+rowOPSmatrix(end-j,:),2);
                    rowOPSmatrix(end-j,:)=mod(rowOPSmatrix(end-j,:)+rowOPSmatrix(r,:),2);
                    rowOPSmatrix(r,:)=mod(rowOPSmatrix(r,:)+rowOPSmatrix(end-j,:),2);
                end
                upperrowsupport=colspan([zeros(1,r),HB.data(r,r+1:end)]>0);%Need to reevaluate upper row support after rowswap
            end
            if ~isempty(upperrowsupport)
                fprintf('Swapping columns %d and %d\n\n',r,upperrowsupport(1));
                HB.data(:,r)=mod(HB.data(:,r)+HB.data(:,upperrowsupport(1)),2);
                HB.data(:,upperrowsupport(1))=mod(HB.data(:,upperrowsupport(1))+HB.data(:,r),2);
                HB.data(:,r)=mod(HB.data(:,r)+HB.data(:,upperrowsupport(1)),2);
                colswaps=[colswaps;[r,upperrowsupport(1)]]; 
                %Keeping track of the column-swaps (via colswaps) of HB to correct kernel
                %basis coordinates so that final kernel corresponds to OG HB
                lowercolsupport=rowspan([zeros(r-1,1);HB.data(r:end,r)]'>0); %Need to reeavuate lowercolsupport after swap
            end
        end
        if length(lowercolsupport)>1 
            %Ellimination of 1s under entry r,r.
            HB.data(lowercolsupport(2:end),:)=mod(HB.data(lowercolsupport(2:end),:)+HB.data(r,:),2);
            %Repeat for rowOPSmatrix
            rowOPSmatrix(lowercolsupport(2:end),:)=mod(rowOPSmatrix(lowercolsupport(2:end),:)+rowOPSmatrix(r,:),2);
        end
    end
    fprintf(strcat('Number of column swaps performed: ',num2str(height(colswaps)),'\n'));
    for c=flip(1:shortside,2) %Upper triangulation
        % c=shortside+1-j;
        if HB.data(c,c)==1
            %Ellimination of 1s above entry c,c.
            uppercolsupport=rowspan([HB.data(1:c,c);zeros(shortside-c,1)]>0);
            if length(uppercolsupport)>1
                HB.data(uppercolsupport(1:end-1),:)=mod(HB.data(uppercolsupport(1:end-1),:)+HB.data(c,:),2);
                %Repeat for rowOPSmatrix
                rowOPSmatrix(uppercolsupport(1:end-1),:)=mod(rowOPSmatrix(uppercolsupport(1:end-1),:)+rowOPSmatrix(c,:),2);
            end
        end
    end
    % Read row-rank of matrix from row-reduced form:
    rowrank=sum(diag(HB.data));
    genmat=ArrayWithHandle();
    genmat.data=[HB.data(1:rowrank,rowrank+1:end)',eye(width(HB.data)-rowrank)];
    %Undo column swaps in the generator matrix.
    for pair=flipud(colswaps)' %more general flip function: flip(A,dim)
        genmat.data(:,pair(1))=mod(genmat.data(:,pair(1))+genmat.data(:,pair(2)),2);
        genmat.data(:,pair(2))=mod(genmat.data(:,pair(2))+genmat.data(:,pair(1)),2);
        genmat.data(:,pair(1))=mod(genmat.data(:,pair(1))+genmat.data(:,pair(2)),2);
        %repeat for colOPSmatrix
        colOPSmatrix(:,pair(1))=mod(colOPSmatrix(:,pair(1))+colOPSmatrix(:,pair(2)),2);
        colOPSmatrix(:,pair(2))=mod(colOPSmatrix(:,pair(2))+colOPSmatrix(:,pair(1)),2);
        colOPSmatrix(:,pair(1))=mod(colOPSmatrix(:,pair(1))+colOPSmatrix(:,pair(2)),2);
    end
end
