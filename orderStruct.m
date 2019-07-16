function s = orderStruct(s1)
  [~, order] = sort(lower(fieldnames(s1)));
  s = orderfields(s1, order);
end
