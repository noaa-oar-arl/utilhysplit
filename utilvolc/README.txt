old clean branch is in /hysplit-users/alicec/clean/utilhysplit. it still has some untracked files in it.

TCM

   volcinverse.py
       classes -
       InvDirectory
       RunInversion
          uses VolcatHysplit paired data class.


MODEL+VOLCAT EVALUATION

   volcpaired.py
       VolcatHysplit class


VOLCAT

    volcat.py for reading VOLCAT netcdf event files.
       DONE - TODO - _get_time function needs checking / updating
       TODO - find_volcat function may need simplifying / splitting.
       TODO - might make some of the functions into a class with an interface.


    volcat_files.py 
       TODO - may be some unused functions in here.

    volcat_plots.py
       TODO - the box plot and cdf plotting functionality needs updating.
              

    volcat_2dplots.py


    TODO - need more functions/methods for identifying/combining volcat files that are at same time
           but have different id numbers. 
           this has been started in volcat_events with the OneEventTime class.

   volcat_event.py

        FUNCTIONS 
            create_event_from_fnames
            create_event 

        Events class

          Attributes:
              self.eventdf - EventDF class

          Properties
              edf - dataframe of self.eventdf 

          Probably has too many methods!

        EventDisplay class (helper class)
        EventStatus  class (helper class)
        EventDF      class (helper class)
        OneEventTime class (helper class)


Model data
   model_event.py contains classes for helping with QVA creation.

        ModelEvent class 
             DataInsertion outputs
             Model runs for input into inversion algorithm.
             Model runs which are result of emissions with inversion algorithm.
             Model runs which are the result of 'traditional' run.

        ModelForecast class
             Instances of the forecast class are created by the ModelEvent class.
             the forecast class should also be able to ingest the VOLCAT data.

   ensemble_polygons contains classes for taking model data and creating polygons

 
qvainterface.py
   TODO - just started writing this.




INTERFACES

    ashapp/ashruninterface  classes for running HYSPLIT model, postprocessing output and creating graphics.

           DONE - the dispersion model runs will run automatically.
           TODO - possibly add some type hinting.
    
    utilvolc/inversioninterface 


DATA INSERTION
    ashapp/DImain.py 
    ashapp/make_volcat_database.py  run this to create a csv file from the sumdf file created by get_summary_file_df.


SOME GUI attempts
    ashapp/volc_ui.py
    ashapp/



