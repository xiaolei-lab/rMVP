#ifndef RMVP_
#define RMVP_
#endif

#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS 1

#include <RcppArmadillo.h>
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <R_ext/Print.h>
#include <progress.hpp>
#include "progress_bar.hpp"
#include "mvp_omp.h"

class MinimalProgressBar_plus: public ProgressBar{
	public:
	MinimalProgressBar_plus()  {
		_finalized = false;
	}

	~MinimalProgressBar_plus() {}

	std::string _time_to_string(double seconds, float progress) {
	  
	    int time = (int) seconds;
	  
	    int hour = 0;
	    int min = 0;
	    int sec = 0;
	  
	    hour = time / 3600;
	    time = time % 3600;
	    min = time / 60;
	    time = time % 60;
	    sec = time;
	  
	    std::stringstream time_strs;
	    time_strs << (progress < 1.0 ? "TimeLeft: " : "RunTime: ");
	    if (hour != 0) time_strs << hour << "h";
	    if (hour != 0 || min != 0) time_strs << min << "m";
	    time_strs << sec << "s";
	    std::string time_str = time_strs.str();
	  
	    return time_str;
	}

	int _compute_nb_ticks(float progress) {
	    return int(progress * _max_ticks);
	}

	std::string _construct_ticks_display_string(int nb) {
	    std::stringstream ticks_strs;
	    for (int i = 1; i <= _max_ticks; ++i) {
	    	if(i < 4){
	    		ticks_strs << ">";
	    	} else if (i < nb) {
	            ticks_strs << "-";
	        } else if(i == nb) {
	        	ticks_strs << ">";
	        } else {
	            ticks_strs << " ";
	        }
	    }
	    std::string tick_space_string = ticks_strs.str();
	    return tick_space_string;
	}

	void end_display() {
      update(1);
    }

    void _finalize_display() {
      if (_finalized) return;
      REprintf("\n");
      _finalized = true;
    }

	void display() {
		flush_console();
	}

	void update(float progress) {
	  
	    // stop if already finalized
	    if (_finalized) return;
	  
	    // start time measurement when update() is called the first time
	    if (_timer_flag) {
	        _timer_flag = false;
	        // measure start time
	        time(&start);
	    } else {
	    
	    	int nb_ticks = _compute_nb_ticks(progress);
		    int delta = nb_ticks - _ticks_displayed;
		    if (delta > 0) {
		    	_ticks_displayed = nb_ticks;
		        std::string cur_display = _construct_ticks_display_string(nb_ticks);
		        
		        // measure current time
		        time(&end);
	    
		        // calculate passed time and remaining time (in seconds)
		        double pas_time = std::difftime(end, start);
		        double rem_time = (progress < 1.0 ? (pas_time / progress) * (1 - progress) : pas_time);
		        if(rem_time < 1 && rem_time > 0.5)	rem_time = 1;

		    	// convert seconds to time string
		        std::string time_string = _time_to_string(rem_time, progress);
		        
		        // ensure overwriting of old time info
		        int empty_length = time_string.length();
		        
	        	std::string empty_space;

		        // merge progress bar and time string
		        std::stringstream strs;
		        if(empty_length_p && abs(empty_length - empty_length_p)){
		        	empty_space = std::string(abs(empty_length - empty_length_p), ' ');
		        	strs << "[" << cur_display << "] " << time_string << empty_space;
		        }else{
		        	strs << "[" << cur_display << "] " << time_string;
		        }
		    	empty_length_p = empty_length;
		    	// strs << "[" << cur_display << "]";

		        std::string temp_str = strs.str();
		        char const* char_type = temp_str.c_str();
		    
		        // print: remove old display and replace with new
		        REprintf("\r");
		        REprintf("%s", char_type);

		    }
		    if (_ticks_displayed >= _max_ticks)
       			// end_display();
       			_finalize_display();
	    }
	}

	void flush_console() {
		#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
			R_FlushConsole();
		#endif
    }

	private: 
	int	empty_length_p = 0;
    int _max_ticks = 45;
    bool _finalized;
    bool _timer_flag = true;
    time_t start, end;
    int _ticks_displayed = 0;
};

class MinimalProgressBar_perc: public ProgressBar{
	public:
	MinimalProgressBar_perc()  {
		_finalized = false;
	}

	~MinimalProgressBar_perc() {}

	std::string _time_to_string(double seconds, float progress) {
	  
	    int time = (int) seconds;
	  
	    int hour = 0;
	    int min = 0;
	    int sec = 0;
	  
	    hour = time / 3600;
	    time = time % 3600;
	    min = time / 60;
	    time = time % 60;
	    sec = time;
	  
	    std::stringstream time_strs;
	    time_strs << (progress < 1.0 ? "TimeLeft: " : "RunTime: ");
	    if (hour != 0) time_strs << hour << "h";
	    if (hour != 0 || min != 0) time_strs << min << "m";
	    time_strs << sec << "s";
	    std::string time_str = time_strs.str();
	  
	    return time_str;
	}

	int _compute_nb_ticks(float progress) {
	    return int(progress * _max_ticks);
	}

	std::string _construct_ticks_display_string(int nb) {
	    std::stringstream ticks_strs;
	    for (int i = 1; i <= _max_ticks; ++i) {
	    	if (i <= nb) {
	            ticks_strs << "*";
	        } else {
	            ticks_strs << " ";
	        }
	    }
	    std::string tick_space_string = ticks_strs.str();
	    return tick_space_string;
	}

	void end_display() {
      update(1);
    }

    void _finalize_display() {
      if (_finalized) return;
      REprintf("\n");
      _finalized = true;
    }

	void display() {
		REprintf("0%%   10   20   30   40   50   60   70   80   90   100%%\n");
      	REprintf("[----|----|----|----|----|----|----|----|----|----|\n");
      	flush_console();
	}

	void update(float progress) {
	  
	    // stop if already finalized
	    if (_finalized) return;
	  
	    // start time measurement when update() is called the first time
	    if (_timer_flag) {
	        _timer_flag = false;
	        // measure start time
	        time(&start);
	    } else {
	    
	    	int nb_ticks = _compute_nb_ticks(progress);
		    int delta = nb_ticks - _ticks_displayed;
		    if (delta > 0) {
		    	_ticks_displayed = nb_ticks;
		        std::string cur_display = _construct_ticks_display_string(nb_ticks);
		        
		        // measure current time
		        time(&end);
	    
		        // calculate passed time and remaining time (in seconds)
		        double pas_time = std::difftime(end, start);
		        double rem_time = (progress < 1.0 ? (pas_time / progress) * (1 - progress) : pas_time);
		        if(rem_time < 1 && rem_time > 0.5)	rem_time = 1;

		    	// convert seconds to time string
		        std::string time_string = _time_to_string(rem_time, progress);
		        
		        // ensure overwriting of old time info
		        int empty_length = time_string.length();
		        
	        	std::string empty_space;

		        // merge progress bar and time string
		        std::stringstream strs;
		        if(empty_length_p && abs(empty_length - empty_length_p)){
		        	empty_space = std::string(abs(empty_length - empty_length_p), ' ');
		        	strs << "[" << cur_display << "] " << time_string << empty_space;
		        }else{
		        	strs << "[" << cur_display << "] " << time_string;
		        }
		    	empty_length_p = empty_length;
		    	// strs << "[" << cur_display << "]";

		        std::string temp_str = strs.str();
		        char const* char_type = temp_str.c_str();
		    
		        // print: remove old display and replace with new
		        REprintf("\r");
		        REprintf("%s", char_type);

		    }
		    if (_ticks_displayed >= _max_ticks)
       			// end_display();
       			_finalize_display();
	    }
	}

	void flush_console() {
		#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
			R_FlushConsole();
		#endif
    }

	private: 
	int	empty_length_p = 0;
    int _max_ticks = 49;
    bool _finalized;
    bool _timer_flag = true;
    time_t start, end;
    int _ticks_displayed = 0;
};

class MinimalProgressBar: public ProgressBar{
	public:
	MinimalProgressBar()  {
		_finalized = false;
	}
	~MinimalProgressBar() {}
	void display() {}
	void update(float progress) {
		if (_finalized) return;
		REprintf("\r");
		REprintf("Calculating in process...(finished %.2f%%)", progress * 100);
	}
	void end_display() {
	if (_finalized) return;
		REprintf("\r");
		
		REprintf("Calculating in process...(finished 100.00%%)");
		REprintf("\n");
		_finalized = true;
	}
	private:
	bool _finalized;
};
