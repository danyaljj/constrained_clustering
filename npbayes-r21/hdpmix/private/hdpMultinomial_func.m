function func = hdpMultinomial_func;

func.iterate   = @hdpMultinomial_iterate;
func.predict   = @hdpMultinomial_predict;
func.numitems  = @Multinomial_numitems;
func.newclass  = @Multinomial_newclass;
func.additem   = @Multinomial_additem;
func.additems  = @Multinomial_additems;
func.delitem   = @Multinomial_delitem;
func.delitems  = @Multinomial_delitems;

