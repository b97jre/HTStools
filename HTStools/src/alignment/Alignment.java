package alignment;
import general.Functions;
import general.RNAfunctions;


public class Alignment {
	
	private int[] start;
	private int[] stop;
	private char[][] Sequence;
	private Alignment nextalignment;
	
	
	Alignment(int []start,int []stop,char[][] Sequence){
		this.start = start;
		this.stop = stop;
		this.Sequence = Sequence;
		nextalignment = null;
		
	}
	
	public void addAlignment(Alignment newAlignment){
		if(nextalignment == null)
			nextalignment = newAlignment;
		else
			nextalignment.addAlignment(newAlignment);
	}
	
	public char[][] getAlignment(int location,int start,int stop, int genome){
		if(isInbetween(location,this.start[genome],this.stop[genome]))
			return getThisAlignment(location,start,stop,genome);
		else if(nextalignment != null)
			return nextalignment.getAlignment(location,start,stop,genome);
		else
			{
			System.out.println("alignment was not found");
			System.out.println(start+" "+stop+" "+location);
			}
			return null;
	}
	
	private char[][]  getThisAlignment(int location,int start,int stop, int genome){
		int alignDir = 1;
		int position = this.start[genome];
		if(this.start[genome] > this.stop[genome]){
			alignDir = -1;
		}

		int locationPointer = 0;
		int leftPosition = start;
		int rightPosition = stop;
		int geneDir = 1;
		if(start > stop){
			leftPosition = stop;
			rightPosition = start;
			geneDir = -1;
		}

		int pointer = 0;
		int count = 0;
		
		while(pointer < Sequence[genome].length && !isInbetween(position + count*alignDir,start,stop)){
			if(Sequence[genome][pointer] == '-')
				pointer++;
			else {
				pointer++;
				count++;
			}
		}
		int leftPos = pointer;
		
		if(pointer == Sequence[genome].length){
			System.out.println("Somtehing wrong with this alignment" );
			System.out.println(this.start[genome]+" "+this.stop[genome] +" "+start+" "+stop);
			return null;
		}
		
		while(pointer < Sequence[genome].length && isInbetween(position + count*alignDir,start,stop)){
			if(Sequence[genome][pointer] == '-')
				pointer++;
			else {
				pointer++;
				count++;
			}
			if(position + count*alignDir == location)
				locationPointer = pointer;
		}
		
		if(locationPointer == 0)
			locationPointer = pointer;
		int rightPos = pointer;
		
		
		boolean easy = true;
		if(locationPointer - leftPos >  ((location-leftPosition)*2) && location-leftPosition > 30){
			easy = false;
		}		
		if(rightPos - locationPointer > ((rightPosition - location)*2)){
			easy = false;
			//System.out.print("Right pos "+ (rightPos - locationPointer)+"\t"+((rightPosition - location)*3));
			//rightPos = locationPointer + ((rightPosition - location)*3);
			//System.out.println("\t"+(rightPos - locationPointer));
		}
		
		//System.out.println(leftPos+ "\t"+rightPos+"\t"+locationPointer);
		
		char[][] newAlignment = null;
		if(easy)
			newAlignment = Functions.getSubmatrix(Sequence,0,Sequence.length,leftPos,rightPos);
		else{
			newAlignment = new char[Sequence.length][0];
			for(int i = 0; i < Sequence.length; i++){
				if(Sequence[i] != null){
					//System.out.println(Sequence[i].length);	
					pointer = locationPointer;
					count = 0;
					while(pointer > 0 && pointer < Sequence[i].length-1 &&count < (location - leftPosition)){
						if(Sequence[i][pointer] != '-')
							count++;
						pointer = pointer - alignDir;
					}
					leftPos = pointer;
					
					count = 0;
					pointer = locationPointer + alignDir;
					while(pointer > 0 && pointer < Sequence[i].length-1 && count < rightPosition - location){
						if(Sequence[i][pointer] != '-')
							count++;
						pointer = pointer + alignDir;
					}
					rightPos = pointer;
					try{
						if(leftPos < rightPos)
							newAlignment[i] = Functions.getSubarray(Sequence[i],leftPos,rightPos);
						else
							newAlignment[i] = Functions.getSubarray(Sequence[i],rightPos,leftPos);
					}
					catch(Exception E){
						E.printStackTrace();
						System.out.println("There was something wrong witht the alignemnt of sequenecs");
						System.out.println("Original length: "+ Sequence[genome].length);
						System.out.println("Sequence length: "+ Sequence[i].length);
						System.out.println("right pos      : "+ rightPos);
						System.out.println("left pos       : "+ leftPos);
						newAlignment[i] = null;
					}
						
				}
				else
					newAlignment[i]=null;
			}
		}
			
		if(geneDir != alignDir)
			newAlignment = RNAfunctions.getComplementary(newAlignment);
		return newAlignment;
	}
	
	private boolean isInbetween(int location, int start, int stop){
		//System.out.println("Start"+start+" stop"+stop);
		if(start < stop){
			if(location >= start && location <=stop)
				return true;
			return false;
		}
		if(location >= stop && location <= start)
			return true;
		return false;
	}
	
	 
}



