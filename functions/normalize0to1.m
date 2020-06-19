function normalized = normalize0to1(mat)
  mat_min = min(mat(:));
  min_max_diff = max(mat(:))-mat_min;
  if min_max_diff == 0 % we do not want to divide by zero
    normalized = zeros(size(mat)); % so return a matrix of zeros because there is no differences in the matrix to normalize
    return
  end
  normalized = (mat-mat_min)./min_max_diff;
end
