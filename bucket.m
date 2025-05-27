function int_rad = bucket(rad,qf,wvl,wvl_rng)

% bin_rng = wvl>wvl_rng(1) & wvl<wvl_rng(2);
% 
% % Integrate over light bucket
% int_rad = zeros(size(rad,2),size(rad,3));
% for i = 1:size(int_rad,1)
%     for j = 1:size(int_rad,2)
% 
%         nv = 0;
%         for n = 1:size(rad,1)
%             if (qf(n,i,j)==0 && bin_rng(n)==1)
%                 int_rad(i,j) = int_rad(i,j) + rad(n,i,j);
%                 nv = nv + 1;
%             end
%         end
% 
%         int_rad(i,j) = int_rad(i,j)*(wvl_rng(2)-wvl_rng(1))/nv;
% 
%     end
% end

bin1 = find(wvl>wvl_rng(1),1);
bin2 = find(wvl<wvl_rng(2),1,'last');

rad(qf~=0) = 0;

int_rad = squeeze(sum(rad(bin1:bin2,:,:),1));
nv = squeeze(sum(qf(bin1:bin2,:,:)==0,1));

int_rad = (wvl_rng(2)-wvl_rng(1))*int_rad./nv;

