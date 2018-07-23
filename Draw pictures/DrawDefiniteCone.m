function DrawDefiniteCone
ezmesh(@(x,z)sqrt(x.*z),[0,1],[0,1])
hold on
ezmesh(@(x,z)-sqrt(x.*z),[0,1],[0,1])
xlabel('x'); ylabel('z'); zlabel('y')
title('y^2=xz');
view([53,26]);
end