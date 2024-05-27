function [translation_vector] = bidirectional_linear(eyePose, params)


%params: direction_azi, direction_ele, factor_azi, factor_ele



trans_azi = params.direction_azi * params.factor_azi * eyePose(1) / 90;
trans_ele = params.direction_ele * params.factor_ele * eyePose(2) / 90;

translation_vector = trans_azi + trans_ele;

end
