#'@keywords internal
definelm <- function(y1,t,df,r,k,seas) {
  p <-  1
  if(r == 1){
    if(k == 1){
      mod <- lm(y1 ~ I(cos(t*(2*pi/seas))^p)   + I(sin(t*(2*pi/seas))^p) +
                  I(t), data = df)
    }else if(k == 2){
      mod <- lm(y1 ~ I(cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                     I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                     I(t  ), data = df)
    }else if(k == 3){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                     I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                     I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                     I(t  ), data = df)
    }else if(k == 4){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                  I(t  ), data = df)
    }else if(k == 5){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p)    + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p)    + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(  cos(t*(2*pi/seas)*4)^p)    + I(  sin(t*(2*pi/seas)*4)^p) +
                  I(  cos(t*(2*pi/seas)*5)^p)    + I(  sin(t*(2*pi/seas)*5)^p) +
                  I(t  ), data = df)
    }else{warning("k is not a positive integer or k > 5")}
  }else if(r == 2){
    if(k == 1){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p) + I( sin(t*(2*pi/seas))^p) +
                     I(t  ) + I((t^2)  ) , data = df)
    }else if(k == 2){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                     I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                     I(t  ) + I((t^2)  ) , data = df)
    }else if(k == 3){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(t  ) + I((t^2)  ) , data = df)
    }else if(k == 4){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                  I(t  ) + I((t^2)  ) , data = df)
    }else if(k == 5){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                  I(  cos(t*(2*pi/seas)*5)^p) + I(  sin(t*(2*pi/seas)*5)^p) +
                  I(t  ) + I((t^2)  ) , data = df)
    }
  }else if(r == 3){
    if(k == 1){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) , data = df)
    }else if(k == 2){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) , data = df)
    }else if(k == 3){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) , data = df)
    }else if(k == 4){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) , data = df)
    }else if(k == 5){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                     I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                     I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                     I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                     I(  cos(t*(2*pi/seas)*5)^p) + I(  sin(t*(2*pi/seas)*5)^p) +
                     I(t  ) + I((t^2)  ) + I((t^3)  ) , data = df)
    }
  }else if(r == 4){
    if(k == 1){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                     I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ), data = df)
    }else if(k == 2){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ), data = df)
    }else if(k == 3){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ), data = df)
    }else if(k == 4){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ), data = df)
    }else if(k == 5){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                     I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                     I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                     I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                     I(  cos(t*(2*pi/seas)*5)^p) + I(  sin(t*(2*pi/seas)*5)^p) +
                     I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ), data = df)
    }
  }else if(r == 5){
    if(k == 1){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ) + I((t^5)  ), data = df)
    }else if(k == 2){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ) + I((t^5)  ), data = df)
    }else if(k == 3){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ) + I((t^5)  ), data = df)
    }else if(k == 4){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                  I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                  I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                  I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                  I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ) + I((t^5)  ), data = df)
    }else if(k == 5){
      mod <- lm(y1 ~ I( cos(t*(2*pi/seas))^p)   + I( sin(t*(2*pi/seas))^p)   +
                     I(  cos(t*(2*pi/seas)*2)^p) + I(  sin(t*(2*pi/seas)*2)^p) +
                     I(  cos(t*(2*pi/seas)*3)^p) + I(  sin(t*(2*pi/seas)*3)^p) +
                     I(  cos(t*(2*pi/seas)*4)^p) + I(  sin(t*(2*pi/seas)*4)^p) +
                     I(  cos(t*(2*pi/seas)*5)^p) + I(  sin(t*(2*pi/seas)*5)^p) +
                     I(t  ) + I((t^2)  ) + I((t^3)  ) + I((t^4)  ) + I((t^5)  ), data = df)
    }
  }else{warning("r is not a positive integer or r > 5")}
  return(mod)
}