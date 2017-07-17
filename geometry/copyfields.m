function copyfields(obj1, obj2)
% copyfields(obj1, obj2)
% Copy the field values of obj2 into the fields in obj1.
% Useful for copying data between handle classes keeping referenced object.

% assert(isequal(class(obj1),class(obj2)),'Different classes')
fields1 = fieldnames(obj1);
fields2 = fieldnames(obj2);
for i = 1 : length(fields1)
  obj1.(fields1{i}) = obj2.(fields2{i});
end

end
