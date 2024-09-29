
for iter = 1:20
    Diff_base_Val = zeros(N, Num_bas_val);
    Int_nu_bu = zeros(N, Num_bas_ue*un);
    Int_cost = zeros(N,1);
    parfor i=1:N
        Diff_base_Val(i,:) = get_basis_Val(Data_save(i+1,1:xn))'-get_basis_Val(Data_save(i,1:xn))';
        int_bu_bu = Data_save(i+1,xn+1:xn+Num_bas_ue^2);
        int_ue_bu = Data_save(i+1,xn+Num_bas_ue^2+1:xn+Num_bas_ue^2+un*Num_bas_ue);
        Qe = Data_save(i+1,xn+Num_bas_ue^2+un*Num_bas_ue+1);
        Ru = int_bu_bu*reshape(Su'*R*Su,[Num_bas_ue^2,1]);
        Int_cost(i,:) = Qe+Ru;
        Int_nu_bu(i,:) = int_ue_bu-int_bu_bu*kron(Su'*R',eye(Num_bas_ue));
    end 
    [ZZ,RankFlag] = chol([Diff_base_Val, 2*Int_nu_bu]'*[Diff_base_Val, 2*Int_nu_bu]);
    RankFlag
    prmt = -pinv([Diff_base_Val, 2*Int_nu_bu])*Int_cost;
    Sv = prmt(1:Num_bas_val);
    Su_pre = Su;
    Su = reshape(prmt(Num_bas_val+1:end),[Num_bas_ue,un])';
end
%save('E:\\OneDrive - nyu.edu\\Code\\ACC_SUMO\\Circular road\\Su.mat','Su')

% NORM = []
% for i = 2:length(SU)
%     NORM = [NORM, norm(SU(:,:,i)-SU(:,:,i-1))];
% end
% plot(1:length(SU)-1, NORM)