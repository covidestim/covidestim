
#Test Case 1- Normal input
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 2 - Negative input
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(-1,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 3 - Date error
date <- as.Date(c("10","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 4 - Character input for Int data
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c("a",0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 5 - Missing Date argument
date <- as.Date(c(""))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 6 - Empty array
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c()
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 7 - Missing Date Argument
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(1,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 8 - Missing integer argument
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 9 - Non dataframe
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- c(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 10 - Check if case is numeric
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- as.integer(c(0,0,0))
hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 11 - Check if hospi is numeric
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
hospi<- as.integer(c(0,0,0))
case <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 12 - Check if deaths is numeric
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
deaths <- as.integer(c(0,0,0))
hospi <- case <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 13 - Check if boost is numeric
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
boost <- as.integer(c(0,0,0))
hospi <- deaths <- case <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 14 - Check if RR is numeric
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
RR<- as.integer(c(0,0,0))
hospi <- deaths <- boost <- case <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 15 - Number of weeks too large
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 7000, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0)

att(validate_df(df))


#Test Case 16 - Test for dataframe interface
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)

cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 


#Test Case 17 - Invalid region
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
deaths <- as.integer(c(0,0,0))
hospi <- case <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Middle-Earth', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 18 - Invalid pop size
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
boost <- as.integer(c(0,0,0))
hospi <- deaths <- case <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Middle-Earth'),start_p_imm=0,cum_p_inf_init=0,df) 

att(validate_df(df))


#Test Case 19 - Invalid start_p_imm
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
RR<- as.integer(c(0,0,0))
hospi <- deaths <- boost <- case <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=-1,cum_p_inf_init=0,df) 


#Test Case 20 - Invalid cum_p_inf_init
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init=-1)


#Test Case 19 - Invalid start_p_imm
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
RR<- as.integer(c(0,0,0))
hospi <- deaths <- boost <- case <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm='a',cum_p_inf_init=0,df) 


#Test Case 20 - Invalid cum_p_inf_init
date <- as.Date(c("2022/10/11","2022/10/12","2022/10/13"))
case <- hospi <- deaths <- boost <- RR <- c(0,0,0)
df <- data.frame(date, case, hospi, deaths, boost, RR)
cfg <- covidestim(nweeks = 31, region = 'Connecticut', pop_size = get_pop('Connecticut'),start_p_imm=0,cum_p_inf_init='a')
