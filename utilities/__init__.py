def is_number(in_value):
    try:
        float(in_value)
        return True
    except ValueError:
        return False
		
def salutation():
    print ("Hello world! This is pyDiffraction library!")