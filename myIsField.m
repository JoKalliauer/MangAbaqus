function isFieldResult = myIsField (inStruct, fieldName)
 %source: https://de.mathworks.com/matlabcentral/answers/95923-is-there-a-matlab-function-that-can-check-if-a-field-exists-in-a-matlab-structure#answer_105275
 % inStruct is the name of the structure or an array of structures to search
 % fieldName is the name of the field for which the function searches
 isFieldResult = false;
 f = fieldnames(inStruct(1));
 for i=1:length(f)
  if(strcmp(f{i},strtrim(fieldName)))
   isFieldResult = true;
   return;
  elseif isstruct(inStruct(1).(f{i}))
   isFieldResult = myIsField(inStruct(1).(f{i}), fieldName);
   if isFieldResult
    return;
   end
  end
end