function [headPoints, pointLabels] = applyEyeTranslation(headPoints, pointLabels, eyePose, translation_method,translation_params)

translation_vector=translation_method(eyePose, translation_params);

headPoints=headPoints+translation_vector;



end