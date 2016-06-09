function xVec = assertVec(x,vecType)

% Created 2/18/16 by DJ (trying to recreate Bryan Conroy's script)

switch vecType
    case 'col'
        xVec = x(:);
    case 'row'
        xVec = x(:)';
end