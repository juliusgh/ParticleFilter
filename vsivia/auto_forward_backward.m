function [x, A19, A20, A21, A22, A23, A24, A25, A26, A27, A28, A29, A30, A31, A32, A33, A34, A35, A36] = auto_forward_backward(x, A19, A20, A21, A22, A23, A24, A25, A26, A27, A28, A29, A30, A31, A32, A33, A34, A35, A36, varargin)

% Forward

B7 = subsref(x, cell2struct({'()' ; { ':'  [1] }}, {'type', 'subs'}, 1)) ;
B8 = subsref(x, cell2struct({'()' ; { ':'  [2] }}, {'type', 'subs'}, 1)) ;
A19 = A19 & (B8+1.0000000000000001e-015) ;
A20 = A20 & (B7/A19) ;
B9 = subsref(x, cell2struct({'()' ; { ':'  [3] }}, {'type', 'subs'}, 1)) ;
B10 = subsref(x, cell2struct({'()' ; { ':'  [4] }}, {'type', 'subs'}, 1)) ;
A21 = A21 & (B10+1.0000000000000001e-015) ;
A22 = A22 & (B9/A21) ;
A23 = A23 & (A20+A22) ;
B11 = subsref(x, cell2struct({'()' ; { ':'  [5] }}, {'type', 'subs'}, 1)) ;
B12 = subsref(x, cell2struct({'()' ; { ':'  [6] }}, {'type', 'subs'}, 1)) ;
A24 = A24 & (B12+1.0000000000000001e-015) ;
A25 = A25 & (B11/A24) ;
A26 = A26 & (A23+A25) ;
B12 = subsref(x, cell2struct({'()' ; { ':'  [6] }}, {'type', 'subs'}, 1)) ;
B7 = subsref(x, cell2struct({'()' ; { ':'  [1] }}, {'type', 'subs'}, 1)) ;
A27 = A27 & (B7+1.0000000000000001e-015) ;
A28 = A28 & (B12/A27) ;
A29 = A29 & (A26+A28) ;
B7 = subsref(x, cell2struct({'()' ; { ':'  [1] }}, {'type', 'subs'}, 1)) ;
B8 = subsref(x, cell2struct({'()' ; { ':'  [2] }}, {'type', 'subs'}, 1)) ;
A30 = A30 & (B7+B8) ;
B9 = subsref(x, cell2struct({'()' ; { ':'  [3] }}, {'type', 'subs'}, 1)) ;
A31 = A31 & (A30+B9) ;
B10 = subsref(x, cell2struct({'()' ; { ':'  [4] }}, {'type', 'subs'}, 1)) ;
A32 = A32 & (A31+B10) ;
B11 = subsref(x, cell2struct({'()' ; { ':'  [5] }}, {'type', 'subs'}, 1)) ;
A33 = A33 & (A32+B11) ;
B12 = subsref(x, cell2struct({'()' ; { ':'  [6] }}, {'type', 'subs'}, 1)) ;
A34 = A34 & (A33+B12) ;
A35 = A35 & (9.9999999999999995e-008*A34) ;
A36 = A36 & (A29+A35) ;

% Backward

A29 = A29 & (A36-A35) ;
A35 = A35 & (A36-A29) ;
A26 = A26 & (A29-A28) ;
A28 = A28 & (A29-A26) ;
A23 = A23 & (A26-A25) ;
A25 = A25 & (A26-A23) ;
A20 = A20 & (A23-A22) ;
A22 = A22 & (A23-A20) ;
B7 = B7 & (A20*A19) ;
A19 = A19 & (B7/A20) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [1] }}, {'type', 'subs'}, 1), B7) ;
B8 = B8 & (A19-1.0000000000000001e-015) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [2] }}, {'type', 'subs'}, 1), B8) ;
B9 = B9 & (A22*A21) ;
A21 = A21 & (B9/A22) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [3] }}, {'type', 'subs'}, 1), B9) ;
B10 = B10 & (A21-1.0000000000000001e-015) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [4] }}, {'type', 'subs'}, 1), B10) ;
B11 = B11 & (A25*A24) ;
A24 = A24 & (B11/A25) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [5] }}, {'type', 'subs'}, 1), B11) ;
B12 = B12 & (A24-1.0000000000000001e-015) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [6] }}, {'type', 'subs'}, 1), B12) ;
B12 = B12 & (A28*A27) ;
A27 = A27 & (B12/A28) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [6] }}, {'type', 'subs'}, 1), B12) ;
B7 = B7 & (A27-1.0000000000000001e-015) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [1] }}, {'type', 'subs'}, 1), B7) ;
A34 = A34 & (A35/9.9999999999999995e-008) ;
A33 = A33 & (A34-B12) ;
B12 = B12 & (A34-A33) ;
A32 = A32 & (A33-B11) ;
B11 = B11 & (A33-A32) ;
A31 = A31 & (A32-B10) ;
B10 = B10 & (A32-A31) ;
A30 = A30 & (A31-B9) ;
B9 = B9 & (A31-A30) ;
B7 = B7 & (A30-B8) ;
B8 = B8 & (A30-B7) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [1] }}, {'type', 'subs'}, 1), B7) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [2] }}, {'type', 'subs'}, 1), B8) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [3] }}, {'type', 'subs'}, 1), B9) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [4] }}, {'type', 'subs'}, 1), B10) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [5] }}, {'type', 'subs'}, 1), B11) ;
x = subsasgn(x, cell2struct({'()' ; { ':'  [6] }}, {'type', 'subs'}, 1), B12) ;

end
