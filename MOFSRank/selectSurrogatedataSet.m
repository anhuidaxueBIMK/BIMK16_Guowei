function [ DB ] = selectSurrogatedataSet( DB,flag )
%SELECTSURROGATEDATASET �˴���ʾ�йش˺�����ժҪ
%   ȥ����ͬԪ��
% dimension = size(surrogateDataSet.TrainX,2);
% surrogateTrainLength = DB.trainSize;
A = [DB.X(1:DB.vernier,:),DB.Y(1:DB.vernier,:) ];
 A = sortrows(A,size(A,2)-1);
[A,remainderSize]=uniqrow(A);
DB.X(1:remainderSize,:) = A(:,1:end-2);
DB.Y(1:remainderSize,:) = A(:,end-1:end);
DB.vernier  = remainderSize;
% if(flag==1)
% %����ģ�����ݿ�
%     if DB.currentSize>DB.trainSize
%         FunctionValue = DB.TrainY(1:DB.currentSize,:);
%         [FrontValue,MaxFront] = F_NDSort(FunctionValue,'all');
%         CrowdDistance         = F_distance(FunctionValue,FrontValue);
%         %
%         N = surrogateTrainLength;
%         Next        = zeros(1,N);
%         while(1)
%             NoN         = numel(FrontValue,FrontValue<MaxFront);%���ط���������Ԫ�ظ���
%             if (NoN>N)
%                 MaxFront = MaxFront-1;
%             else
%                 break;
%             end
%         end
%        Next(1:NoN) = find(FrontValue<MaxFront);
%         %
% %         Last          = find(FrontValue==MaxFront);
% %         [~,Rank]      = sort(CrowdDistance(Last),'descend');
% %         Next(NoN+1:N) = Last(Rank(1:N-NoN));
%     DB.TrainX(1:NoN,:) = DB.TrainX(Next(1:NoN),:);
%     DB.TrainY(1:NoN,:) = DB.TrainY(Next(1:NoN),:);
%     DB.currentSize  = NoN;
%     end
% end
end

