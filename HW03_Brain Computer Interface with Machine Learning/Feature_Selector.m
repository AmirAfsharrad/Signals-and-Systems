function [TrainSet, labels, feature_list, TestSet] = Feature_Selector(FData, exe_or_img, number_of_features)
    if (strcmp(exe_or_img, 'exe'))
        [J(:,1),~,Feature1] = Jvalue(FData.exe.Arm, FData.exe.Leg);
        [J(:,2),~] = Jvalue(FData.exe.Arm, FData.exe.Thumb);
        [J(:,3),~] = Jvalue(FData.exe.Arm, FData.exe.Idle);
        [J(:,4),~] = Jvalue(FData.exe.Leg, FData.exe.Thumb);
        [J(:,5),~] = Jvalue(FData.exe.Leg, FData.exe.Idle);
        [J(:,6),J_Index,Feature2] = Jvalue(FData.exe.Thumb, FData.exe.Idle);
        
        [~,~,TestFeature] = Jvalue(FData.exe.test, FData.exe.test);
        n = size(TestFeature,1);
        TestFeature = TestFeature(1:n/2,:);
        
        Feature = [Feature1; Feature2];
        clear Feature1 Feature2
        
        [~,I] = sort(J,1,'descend');
        
        N = floor(number_of_features/6)*ones(1,6);
        remainder = mod(number_of_features,6);
        N(1:remainder) = N(1:remainder)+1;
        
        TestSet(:,1:N(1)) = TestFeature(:,I(1:N(1),1));
        TestSet(:,sum(N(1:1))+1:sum(N(1:2))) = TestFeature(:,I(1:N(2),2));
        TestSet(:,sum(N(1:2))+1:sum(N(1:3))) = TestFeature(:,I(1:N(3),3));
        TestSet(:,sum(N(1:3))+1:sum(N(1:4))) = TestFeature(:,I(1:N(4),4));
        TestSet(:,sum(N(1:4))+1:sum(N(1:5))) = TestFeature(:,I(1:N(5),5));
        TestSet(:,sum(N(1:5))+1:sum(N(1:6))) = TestFeature(:,I(1:N(6),6));
        
        
        TrainSet(:,1:N(1)) = Feature(:,I(1:N(1),1));
        TrainSet(:,sum(N(1:1))+1:sum(N(1:2))) = Feature(:,I(1:N(2),2));
        TrainSet(:,sum(N(1:2))+1:sum(N(1:3))) = Feature(:,I(1:N(3),3));
        TrainSet(:,sum(N(1:3))+1:sum(N(1:4))) = Feature(:,I(1:N(4),4));
        TrainSet(:,sum(N(1:4))+1:sum(N(1:5))) = Feature(:,I(1:N(5),5));
        TrainSet(:,sum(N(1:5))+1:sum(N(1:6))) = Feature(:,I(1:N(6),6));
        
        labels(1:20) = 1;
        labels(21:40) = 2;
        labels(41:60) = 3;
        labels(61:80) = 4;
        
        feature_list(1:N(1),:) = J_Index(I(1:N(1),1),:);
        feature_list(sum(N(1:1))+1:sum(N(1:2)),:) = J_Index(I(1:N(2),2),:);
        feature_list(sum(N(1:2))+1:sum(N(1:3)),:) = J_Index(I(1:N(3),3),:);
        feature_list(sum(N(1:3))+1:sum(N(1:4)),:) = J_Index(I(1:N(4),4),:);
        feature_list(sum(N(1:4))+1:sum(N(1:5)),:) = J_Index(I(1:N(5),5),:);
        feature_list(sum(N(1:5))+1:sum(N(1:6)),:) = J_Index(I(1:N(6),6),:);
    end
    
	if (strcmp(exe_or_img, 'img'))
        [J(:,1),~,Feature1] = Jvalue(FData.img.Arm, FData.img.Leg);
        [J(:,2),~] = Jvalue(FData.img.Arm, FData.img.Thumb);
        [J(:,3),J_Index,Feature2] = Jvalue(FData.img.Leg, FData.img.Thumb);
        
        [~,~,TestFeature] = Jvalue(FData.img.test, FData.img.test);
        n = size(TestFeature,1);
        TestFeature = TestFeature(1:n/2,:);
        
        Feature = [Feature1; Feature2(21:40,:)];
        clear Feature1 Feature2
        
        [~,I] = sort(J,1,'descend');
        
        N = floor(number_of_features/3)*ones(1,3);
        remainder = mod(number_of_features,3);
        N(1:remainder) = N(1:remainder)+1;
        
        TestSet(:,1:N(1)) = TestFeature (:,I(1:N(1),1));
        TestSet(:,sum(N(1:1))+1:sum(N(1:2))) = TestFeature (:,I(1:N(2),2));
        TestSet(:,sum(N(1:2))+1:sum(N(1:3))) = TestFeature (:,I(1:N(3),3));
        
        TrainSet(:,1:N(1)) = Feature(:,I(1:N(1),1));
        TrainSet(:,sum(N(1:1))+1:sum(N(1:2))) = Feature(:,I(1:N(2),2));
        TrainSet(:,sum(N(1:2))+1:sum(N(1:3))) = Feature(:,I(1:N(3),3));
        
        labels(1:20) = 1;
        labels(21:40) = 2;
        labels(41:60) = 3;
        
        feature_list(1:N(1),:) = J_Index(I(1:N(1),1),:);
        feature_list(sum(N(1:1))+1:sum(N(1:2)),:) = J_Index(I(1:N(2),2),:);
        feature_list(sum(N(1:2))+1:sum(N(1:3)),:) = J_Index(I(1:N(3),3),:);
	end

end