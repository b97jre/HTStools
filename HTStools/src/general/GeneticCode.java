package general;

public class GeneticCode {

	private char[][][] RNA2protein;


	public GeneticCode(){
		RNA2protein = new char [6][6][6];
		for(int i= 0; i < 6;i++){
			for(int j= 0; j < 6;j++){
				for(int k= 0; k < 6;k++){
					RNA2protein[i][j][k]='n';
				}
			}
		}
		RNA2protein[1][1][1]='K';
		RNA2protein[1][1][2]='N';
		RNA2protein[1][1][3]='K';
		RNA2protein[1][1][4]='N';
		RNA2protein[1][2][1]='T';
		RNA2protein[1][2][2]='T';
		RNA2protein[1][2][3]='T';
		RNA2protein[1][2][4]='T';
		RNA2protein[1][3][1]='R';
		RNA2protein[1][3][2]='S';
		RNA2protein[1][3][3]='R';
		RNA2protein[1][3][4]='S';
		RNA2protein[1][4][1]='I';
		RNA2protein[1][4][2]='I';
		RNA2protein[1][4][3]='M';
		RNA2protein[1][4][4]='I';
		RNA2protein[2][1][1]='Q';
		RNA2protein[2][1][2]='H';
		RNA2protein[2][1][3]='Q';
		RNA2protein[2][1][4]='H';
		RNA2protein[2][2][1]='P';
		RNA2protein[2][2][2]='P';
		RNA2protein[2][2][3]='P';
		RNA2protein[2][2][4]='P';
		RNA2protein[2][3][1]='R';
		RNA2protein[2][3][2]='R';
		RNA2protein[2][3][3]='R';
		RNA2protein[2][3][4]='R';
		RNA2protein[2][4][1]='L';
		RNA2protein[2][4][2]='L';
		RNA2protein[2][4][3]='L';
		RNA2protein[2][4][4]='L';
		RNA2protein[3][1][1]='E';
		RNA2protein[3][1][2]='D';
		RNA2protein[3][1][3]='E';
		RNA2protein[3][1][4]='D';
		RNA2protein[3][2][1]='A';
		RNA2protein[3][2][2]='A';
		RNA2protein[3][2][3]='A';
		RNA2protein[3][2][4]='A';
		RNA2protein[3][3][1]='G';
		RNA2protein[3][3][2]='G';
		RNA2protein[3][3][3]='G';
		RNA2protein[3][3][4]='G';
		RNA2protein[3][4][1]='V';
		RNA2protein[3][4][2]='V';
		RNA2protein[3][4][3]='V';
		RNA2protein[3][4][4]='V';
		RNA2protein[4][1][1]='*';
		RNA2protein[4][1][2]='Y';
		RNA2protein[4][1][3]='*';
		RNA2protein[4][1][4]='Y';
		RNA2protein[4][2][1]='S';
		RNA2protein[4][2][2]='S';
		RNA2protein[4][2][3]='S';
		RNA2protein[4][2][4]='S';
		RNA2protein[4][3][1]='*';
		RNA2protein[4][3][2]='C';
		RNA2protein[4][3][3]='W';
		RNA2protein[4][3][4]='C';
		RNA2protein[4][4][1]='L';
		RNA2protein[4][4][2]='F';
		RNA2protein[4][4][3]='L';
		RNA2protein[4][4][4]='F';	
	}
	
	public char TranslateTriplet(int A,int B,int C){
		return RNA2protein[A][B][C];
	}

	public char[] TranslateRNAseq(int start, int [] sequence){
		int pointer =  start;
		char []proteinSeq = new char[sequence.length/3];
		int proteinLength = 0;
		
		while(sequence.length > (pointer+2) && RNA2protein[sequence[pointer]][sequence[pointer+1]][sequence[pointer+2]]!='*'){
			proteinSeq[proteinLength]=RNA2protein[sequence[pointer]][sequence[pointer+1]][sequence[pointer+2]];
//			System.out.print(sequence[pointer]+" "+sequence[pointer+1]+" "+sequence[pointer+2]);
//			System.out.println("="+proteinSeq[proteinLength]);
			proteinLength++;
			pointer=pointer+3;
		}
		if(sequence.length > pointer+2){
			proteinSeq[proteinLength]='*';
			proteinLength++;
			pointer=pointer+3;
		}
		char [] FS= new char[proteinLength];
		for(int i = 0; i < proteinLength;i++){
			FS[i] = proteinSeq[i];
		}
		return FS;
		
		
		
	}
}




