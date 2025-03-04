function [translation_vector] = bidirectional_declining_sin(eyePose, params)

% params: direction_azi,direction_ele,factor_azi,factor_ele

trans_azi = params.direction_azi * params.factor_azi * sin(deg2rad(eyePose(1)));
trans_ele = params.direction_ele * params.factor_ele * sin(deg2rad(eyePose(2)));

translation_vector = trans_azi + trans_ele;
    
end
    