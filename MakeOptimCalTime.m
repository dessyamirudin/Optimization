%Model Optimisation under Uncertainty

%Author  Dessy Amirudin
%Master Student
%Operational Research with Computational Optimization
%The University of Edinburgh

%terminating parameter
epsilon=1e-6;
P= [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 0 0 0];

Q=inv(eye(4)-P');

%C1 is buffer 1 and 3
%C2 is buffer 2 and 4
C=[1 3;
    2 4];
[numclass,numbuffclass]=size(C);
numbuff=numclass*numbuffclass;

%STATISTIK
numberrealization=[];
minOptim=[];
timeOptim=[];
combinationreal=[];
optimVkreal=[];

%number of realization
maxNumReal=2;

for NumReal=2:maxNumReal
%combination of realization
matrcomb=dec2bin(0:2^NumReal-1,2)-'0';
posreal=matrcomb;
[rowC,colC]=size(posreal);

for n=1:rowC
    for m=1:colC
        posreal(n,m)=posreal(n,m)+1;
    end
end

a=[30 15 20 10];
m=[0.1 0.2 0.05 0.3];
%m=0.4*m;
tr=10;
tr=7*tr;


%probability of each realization, the value of probability is same for each
%realization
probreal=[];
for i=1:NumReal
    probreal=[probreal 1/NumReal];
end

%alpha on each realization
lbalpha=0.5;
ubalpha=1.5;
alpha=[];

for i=0:NumReal-1
    alpha=[alpha;i*(ubalpha-lbalpha)/(NumReal-1)+lbalpha zeros(1,3)];
end

%using simplex method
options = optimset('LargeScale', 'off', 'Simplex', 'on');

%matrix for Optimum Value. I use the upper bound as basis of optimal value
OptimValue=[];

%Matrix to VK Optimum
Combination=[];
OptimVK=[];

%start time calculation
tStart=tic;

%$$$$$$$$$$$$$$ Calculating Boundary $$$$$$$$$$$$$$

A1=[];
B1=[];
A2=[];
B2=[];
f=[];
boundcalc=[];
boundlimit=1.2;

for indexreal=1:NumReal
    A1nom=[];
    B1den=[];
    for n=1:4
        A1nom=[A1nom;Q(n,:)*transpose(a)*m(n)];
        B1den=[B1den;Q(n,:)*transpose(alpha(indexreal,:))*m(n)];
    end
    A1=[A1 A1nom];
    B1=[B1 B1den];
    
    A2nom=[];
    for i=1:2
        sumA=0;
        for j=1:2
            sumA=sumA+Q(C(i,j),:)*(a+tr*alpha(indexreal,:))'*m(C(i,j));
        end
        A2nom=[A2nom;sumA];
    end
    
    A2=[A2 A2nom];
end

v1max=[];
v1max=max(B1,[],2);

Bden1=[v1max(1)*tr+v1max(3)*tr;
    v1max(2)*tr+v1max(4)*tr];

for border=1:NumReal
    boundcalcx=[];
    for baris=1:2
        boundcalcx=[boundcalcx A2(baris,border)/tr];
    end
    
    boundcalc=[boundcalc;boundcalcx];
end

boundmax=max(boundcalc,[],2);
calccom=0;

for boundcom=1:NumReal
    if boundmax(boundcom)<=boundlimit
        calccom=calccom+1;
    end
end

boundcomposition=calccom/NumReal;

dif=0;
dif1=boundlimit-min(boundmax);
dif2=max(boundmax)-boundlimit;

if dif2>dif1
    dif=1;
end
%$$$$$$$$$$$$$$$$ end boundary calculation

for com=1:rowC
    %reset the diff, LB and UB
    LB=0;
    UB=10000;
    diff=UB-LB;
    realizationonecalc=0;
    
    for indexcom=1:colC
        if posreal(com,indexcom)==1
            realizationonecalc=realizationonecalc+1;
        end
    end
    
    realizationcomposition=realizationonecalc/NumReal;
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Equation One Dominance $$$$$$$$$$$$$$$$$$$$$$
    if realizationcomposition>=0.5 && boundcomposition>=0.5 && dif==0
        fprintf('*****************************************************\n')
        fprintf('*************** COMBINATION :')
        for indexcom=1:colC
            fprintf(' %d',posreal(com,indexcom))
        end
        fprintf(' *****************\n')
        fprintf('*****************************************************\n')
        fprintf('\n')
        Combination=[Combination;posreal(com,:)];
        
        %reset all matrix for first second and third realization
        %matrix for first realization
        fA1=[];
        fB1=[];
        fA2=[];
        fB2=[];
        f=[];
        grad=zeros(numbuff,NumReal);
        addSecond=[];
        
        %######################### ALL REALIZATION ###########################
        
        %constructing vector of nominator and denominator for all realzation
        
        for indexreal=1:NumReal
            if posreal(com,indexreal)==1
                fA2=[fA2 zeros(numclass,1)];
                fB2=[fB2 zeros(numclass,1)];
                addSecond=[addSecond;zeros(2,4)];
                %calculating nominator and denominator first realization
                %denominator is equal to the lower bound of vk
                fA1nom=[];
                fB1den=[];
                for n=1:4
                    fA1nom=[fA1nom;Q(n,:)*transpose(a)*m(n)];
                    fB1den=[fB1den;Q(n,:)*transpose(alpha(indexreal,:))*m(n)];
                end
                
                fA1=[fA1 fA1nom];
                fB1=[fB1 fB1den];
            else
                fA1=[fA1 zeros(numbuff,1)];
                fB1=[fB1 zeros(numbuff,1)];
                
                fA2nom=[];
                fB2den=[];
                for i=1:2
                    sumA=0;
                    sumB=0;
                    for j=1:2
                        sumA=sumA+Q(C(i,j),:)*(a+tr*alpha(indexreal,:))'*m(C(i,j));
                        sumB=sumB+Q(C(i,j),:)*alpha(indexreal,:)'*m(C(i,j));
                    end
                    %update matrix
                    fA2nom=[fA2nom;sumA/(1-sumB)];
                    fB2den=[fB2den;sumB];
                    
                end
                
                fA2=[fA2 fA2nom];
                fB2=[fB2 fB2den];
                %calculating the vector that need to be added to first equation
                addFirst=zeros(2,4);
                for i=1:2
                    for j=1:2
                        addFirst(i,C(i,j))=m(C(i,j))/(1-fB2(i,indexreal));
                    end
                end
                
                addSecond=[addSecond;addFirst];
                %masih problem disini, bagaimana cara mengakses sebagian dari
                %isi matriz? misal dari 10 x 10 saya hanya ingin mengakses line
                %8 dan 9
                %%%Solved by (a:b,)
                
            end
        end
        
        %calculating the maximum value of each vk
        vmax=[];
        vmax=max(fB1,[],2)+0.01;
        
        %matrix for additional lower bound
        Alb4=[];
        blb4=[];
        
        %component on the lower bound calculation that fixed
        costlb=[1;0;0;0;0];
        Alb1=[0 1 0 1 0;
            0 0 1 0 1];
        Alb2=[zeros(4,1) -eye(4)];
        Alb3=[-1 zeros(1,4)];
        blb1=ones(2,1);
        blb2=-vmax;
        blb3=0;
        lblb=zeros(5,1);
        Alb=[Alb1;Alb2;Alb3;Alb4];
        blb=[blb1;blb2;blb3;blb4];
        
        count=1;
        
        while diff>epsilon
            %fprintf('LOOP:%d \n',count)
            
            %calculating Lower Bound of the equation
            Alb=[Alb;Alb4];
            blb=[blb;blb4];
            Alb4=[];
            blb4=[];
            %simplex method for Lower Bound
            [varlb,thetalb,exitflaglb,outputlb,Pilb]= linprog(costlb,Alb,blb,[],[],lblb,[],[],options);
            %updating value of LowBound
            if thetalb>=LB
                LB=thetalb;
            end
            %fprintf('Lower Bound :%d \n',thetalb)
            
            %updating the value of vkmax
            vmaxupt=[];
            for n=1:4
                vmaxupt=[vmaxupt;varlb(n+1)];
            end
            
            %reset matrix thetastar...because need to be updated/changed in
            %each iteration
            thetastar=[];
            
            
            for index=1:NumReal
                %********** IF EQUATION 1 REALIZED ****************
                %calculating the value of first,second and third realization with vkmax
                %if the realization chose equation number 1
                if posreal(com,index)==1
                    f=[];
                    for n=1:4
                        f=[f;fA1(n,index)/(vmaxupt(n)-fB1(n,index))];
                    end
                    %finding which (index) vk giving highest value. also
                    %calculating the thetastar
                    [thetastar(index) ind]=max(f);
                    
                    %gradien for first realization
                    grad1=zeros(4,1);
                    grad1(ind)=-fA1(ind,index)/(vmaxupt(ind)-fB1(ind,index))^2;
                    
                    %updating matrix gradien
                    grad(:,index)=grad1;
                else
                    %********** IF EQUATION 2 REALIZED ****************
                    %calculating theta for second relization
                    costvar1real=[1;0;0;0;0];
                    A11=[zeros(4,1) diag(m)];
                    A12=[zeros(4,1) eye(4)-transpose(P)];
                    A13=[ones(2,1) addSecond(index*2-1:index*2,:)];
                    addSecond(index*2-1:index*2,:);
                    %baru sampai
                    %sini-----------------------------
                    b11=vmaxupt*tr;
                    b12=a'+alpha(index,:)'*tr;
                    b13=fA2(:,index)+ones(2,1)*tr;
                    A1=[A11;A12;-A13];
                    b1=[b11;b12;-b13];
                    %simplex method for theta2
                    [var2,theta2,exitflag2,output2,Pi2]= linprog(costvar1real,A1,b1,[],[],[],[],[],options);
                    Pi2.ineqlin;
                    thetastar(index)=theta2;
                    
                    %gradien for second realization
                    grad2=zeros(4,1);
                    for n=1:4
                        grad2(n)=Pi2.ineqlin(n)*tr;
                    end
                    grad(:,index)=grad2;
                    
                end
            end
            thetastar;
            grad';
            
            %************ FINISHING UP and ADDING CUT ******************
            %summing all gradien
            sumgrad2=probreal*grad';
            
            UpB=probreal*thetastar';
            if UpB<UB
                UB=UpB;
            end
            %calculating Upper Bound
            %fprintf('UpB :%d \n',UpB)
            %fprintf('Optimum Value :%d \n',UB)
            diff=UB-LB;
            %diff=0;
            %fprintf('Difference :%d \n',diff)
            %adding new cut
            Alb4=[-1 sumgrad2];
            blb4=-UpB+sumgrad2*vmaxupt;
            
            count=count+1;
            fprintf('\n')
        end
        %Save the optimum value
        OptimValue=[OptimValue;UB];
        %Save the vk for each Optimum Value
        OptimVK=[OptimVK;vmaxupt'];
    end
    
    %$$$$$$$$$$$$$$$$$$$$$$$ Equation Two Dominance $$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    if realizationcomposition<=0.5 && boundcomposition<=0.9 && dif==1
        fprintf('*****************************************************\n')
        fprintf('*************** COMBINATION :')
        for indexcom=1:colC
            fprintf(' %d',posreal(com,indexcom))
        end
        fprintf(' *****************\n')
        fprintf('*****************************************************\n')
        fprintf('\n')
                    
        Combination=[Combination;posreal(com,:)];
        
        %reset all matrix for first second and third realization
        %matrix for first realization
        fA1=[];
        fB1=[];
        fA2=[];
        fB2=[];
        f=[];
        grad=zeros(numbuff,NumReal);
        addSecond=[];
        
        %######################### ALL REALIZATION ###########################
        
        %constructing vector of nominator and denominator for all realzation
        
        for indexreal=1:NumReal
            if posreal(com,indexreal)==1
                fA2=[fA2 zeros(numclass,1)];
                fB2=[fB2 zeros(numclass,1)];
                addSecond=[addSecond;zeros(2,4)];
                %calculating nominator and denominator first realization
                %denominator is equal to the lower bound of vk
                fA1nom=[];
                fB1den=[];
                for n=1:4
                    fA1nom=[fA1nom;Q(n,:)*transpose(a)*m(n)];
                    fB1den=[fB1den;Q(n,:)*transpose(alpha(indexreal,:))*m(n)];
                end
                
                fA1=[fA1 fA1nom];
                fB1=[fB1 fB1den];
            else
                fA1=[fA1 zeros(numbuff,1)];
                fB1=[fB1 zeros(numbuff,1)];
                
                fA2nom=[];
                fB2den=[];
                for i=1:2
                    sumA=0;
                    sumB=0;
                    for j=1:2
                        sumA=sumA+Q(C(i,j),:)*(a+tr*alpha(indexreal,:))'*m(C(i,j));
                        sumB=sumB+Q(C(i,j),:)*alpha(indexreal,:)'*m(C(i,j));
                    end
                    %update matrix
                    fA2nom=[fA2nom;sumA/(1-sumB)];
                    fB2den=[fB2den;sumB];
                    
                end
                
                fA2=[fA2 fA2nom];
                fB2=[fB2 fB2den];
                %calculating the vector that need to be added to first equation
                addFirst=zeros(2,4);
                for i=1:2
                    for j=1:2
                        addFirst(i,C(i,j))=m(C(i,j))/(1-fB2(i,indexreal));
                    end
                end
                
                addSecond=[addSecond;addFirst];
                %masih problem disini, bagaimana cara mengakses sebagian dari
                %isi matriz? misal dari 10 x 10 saya hanya ingin mengakses line
                %8 dan 9
                %%%Solved by (a:b,)
                
            end
        end
        
        %calculating the maximum value of each vk
        vmax=[];
        vmax=max(fB1,[],2)+0.01;
        
        %matrix for additional lower bound
        Alb4=[];
        blb4=[];
        
        %component on the lower bound calculation that fixed
        costlb=[1;0;0;0;0];
        Alb1=[0 1 0 1 0;
            0 0 1 0 1];
        Alb2=[zeros(4,1) -eye(4)];
        Alb3=[-1 zeros(1,4)];
        blb1=ones(2,1);
        blb2=-vmax;
        blb3=0;
        lblb=zeros(5,1);
        Alb=[Alb1;Alb2;Alb3;Alb4];
        blb=[blb1;blb2;blb3;blb4];
        
        count=1;
        
        while diff>epsilon
            %fprintf('LOOP:%d \n',count)
            
            %calculating Lower Bound of the equation
            Alb=[Alb;Alb4];
            blb=[blb;blb4];
            Alb4=[];
            blb4=[];
            %simplex method for Lower Bound
            [varlb,thetalb,exitflaglb,outputlb,Pilb]= linprog(costlb,Alb,blb,[],[],lblb,[],[],options);
            %updating value of LowBound
            if thetalb>=LB
                LB=thetalb;
            end
            %fprintf('Lower Bound :%d \n',thetalb)
            
            %updating the value of vkmax
            vmaxupt=[];
            for n=1:4
                vmaxupt=[vmaxupt;varlb(n+1)];
            end
            
            %reset matrix thetastar...because need to be updated/changed in
            %each iteration
            thetastar=[];
            
            
            for index=1:NumReal
                %********** IF EQUATION 1 REALIZED ****************
                %calculating the value of first,second and third realization with vkmax
                %if the realization chose equation number 1
                if posreal(com,index)==1
                    f=[];
                    for n=1:4
                        f=[f;fA1(n,index)/(vmaxupt(n)-fB1(n,index))];
                    end
                    %finding which (index) vk giving highest value. also
                    %calculating the thetastar
                    [thetastar(index) ind]=max(f);
                    
                    %gradien for first realization
                    grad1=zeros(4,1);
                    grad1(ind)=-fA1(ind,index)/(vmaxupt(ind)-fB1(ind,index))^2;
                    
                    %updating matrix gradien
                    grad(:,index)=grad1;
                else
                    %********** IF EQUATION 2 REALIZED ****************
                    %calculating theta for second relization
                    costvar1real=[1;0;0;0;0];
                    A11=[zeros(4,1) diag(m)];
                    A12=[zeros(4,1) eye(4)-transpose(P)];
                    A13=[ones(2,1) addSecond(index*2-1:index*2,:)];
                    addSecond(index*2-1:index*2,:);
                    %baru sampai
                    %sini-----------------------------
                    b11=vmaxupt*tr;
                    b12=a'+alpha(index,:)'*tr;
                    b13=fA2(:,index)+ones(2,1)*tr;
                    A1=[A11;A12;-A13];
                    b1=[b11;b12;-b13];
                    %simplex method for theta2
                    [var2,theta2,exitflag2,output2,Pi2]= linprog(costvar1real,A1,b1,[],[],[],[],[],options);
                    Pi2.ineqlin;
                    thetastar(index)=theta2;
                    
                    %gradien for second realization
                    grad2=zeros(4,1);
                    for n=1:4
                        grad2(n)=Pi2.ineqlin(n)*tr;
                    end
                    grad(:,index)=grad2;
                    
                end
            end
            thetastar;
            grad';
            
            %************ FINISHING UP and ADDING CUT ******************
            %summing all gradien
            sumgrad2=probreal*grad';
            
            UpB=probreal*thetastar';
            if UpB<UB
                UB=UpB;
            end
            %calculating Upper Bound
            %fprintf('UpB :%d \n',UpB)
            %fprintf('Optimum Value :%d \n',UB)
            diff=UB-LB;
            %diff=0;
            %fprintf('Difference :%d \n',diff)
            %adding new cut
            Alb4=[-1 sumgrad2];
            blb4=-UpB+sumgrad2*vmaxupt;
            
            count=count+1;
            fprintf('\n')
        end
        %Save the optimum value
        OptimValue=[OptimValue;UB];
        %Save the vk for each Optimum Value
        OptimVK=[OptimVK;vmaxupt'];
    end
end
tElapse=toc(tStart);
fprintf('\n')
fprintf('Algorithm Length Calculation: %d',tElapse)
fprintf(' seconds\n')
fprintf('\n')
fprintf('SUMMARY:%d \n')
%Combination=posreal;
Combination;
OptimValue;
[minOptim(NumReal) indOptim(NumReal)]=min(OptimValue);
optimVkreal=[optimVkreal;OptimVK(indOptim(NumReal),:)];
combinationreal=[combinationreal;zeros(1,maxNumReal-NumReal) Combination(indOptim(NumReal),:)];
numberrealization=[numberrealization;NumReal];
timeOptim=[timeOptim;tElapse];
end

numberrealization
minOptim
timeOptim
combinationreal
optimVkreal