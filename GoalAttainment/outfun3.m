function stop = outfun3(x,optimValues,state,history)
    stop = false;
    history.add(x,optimValues.fval);
end
