function s = orderStruct(s1)
%% s = orderStruct(s1)
%
%  Function to order fields in structure alphabetically

[~, order] = sort(lower(fieldnames(s1)));
s = orderfields(s1, order);
end
