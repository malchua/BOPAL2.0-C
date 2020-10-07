
//cette classe permet de faire des alignements de couple de s�quences et tenant compte de seulement les duplications et pertes


//Remarque : dans la procedure de construction de l'ancetre,

//si on arrive a une inversion ou une substitution dans XY, on doit parcourir Y jusqu'a cette inversion



import java.util.ArrayList;
import java.util.Scanner;
import java.io.*;


public class Aligning {


	static String voisin="c";
	
	static String [] s4;
	
	static String [] s5;
	
	static String [] v;
	
	static String ancetre="";
	
	static int cost; // cout du dernier alignement
	
	static boolean show_time=false;
	static boolean show_align=false;
	
	static long time=0;
	static String hhmmss; // temps de calcul sous la forme hh:mm:ss
	
	
	//______________________________________________________________________________
	//|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|
	//|:=--=:|     Main    |:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|:=--=:|

	//Preliminaries: __    __    __    __    __    __    __    __    __    __    __
	//__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \_
	
	
	public static void main(String[] args) {
		
		
		String usage = "Usage: Aligning  -dt  genome1  genome2  neighbor[Optional]"+"\n"+
		"-d show detailed results."+"\n"+
		"-t show running time."+"\n"+
		"genome1, genome2 and neighbor are strings of genes."+"\n"+
		"gene names are strings (delimeted by ',')."+"\n\n"+
		"Example : Java OrthoAlign  -d  abc,dfe,abc  abc,fdsg,cvd  abc,dfe,abc,cvd ."+"\n";
		
		
		if ((args.length<2) || (args.length>4)) System.out.println(usage);
		
		else 
		
		{
			String opt=args[0];
			
			if (opt.equals("-dt") || opt.equals("-d") || opt.equals("-t"))
			
				{
					if (args.length<3) 
						
						System.out.println(usage);
					
					else 
					
					{
						
						if (opt.equals("-dt"))
							
						{
							show_time=true;
							show_align=true;
							
						}
						
						else if (opt.equals("-d"))
							
							show_align=true;
						
						else if (opt.equals("-t"))
							
							show_time=true;
						
						
						
						if (args.length==3) prog1_perte(args[1],args[2],true,true,false,true);
						
						else 
							
						{
							voisin=args[3];
							
							prog1_perte(args[1],args[2],true,false,true,true);
							
						}
						
					}
				
					
				}
				
			
			else 
				
					{
				
						if (args.length==2) prog1_perte(args[0],args[1],true,true,false,true);
						
						else 
							
							
						{
							voisin=args[2];
							
							prog1_perte(args[0],args[1],true,false,true,true);
				
						}
						
				
					}
			
				
		}
		
		
		
		
	}
		
	
	

public static void prog1_perte(String s2, String s3, boolean ps, boolean dp,boolean vois, boolean rp)

{
	
	// meme chose que prog1 sauf qu'il prend en compte des perte --> dup ( al.post_traitement_voisin_perte(C); )
	

s4=s2.split(", *");

s5=s3.split(", *");

if (vois) v=voisin.split(", *");

else v="pas, de. voisin".split(", *");

ancetre="";

cost=0;  // reinitialiser le cout de l'alignement � 0


int []X=new int[s4.length];
int []Y=new int[s5.length];

int []C=new int[v.length];



for (int i = 0; i < X.length; i++) 
	
	{
	    if (s4[i].charAt(0)=='-')
	    
		X[i] = -(s4[i].substring(1).hashCode());
		
		else X[i] = (s4[i].hashCode());
		
		
		
	}



for (int i = 0; i < Y.length; i++) 
	
{
	
	if (s5[i].charAt(0)=='-')
	
		Y[i] = -(s5[i].substring(1).hashCode());	
		
	else Y[i] = (s5[i].hashCode());
	
	
	
}
	


for (int i = 0; i < C.length; i++) 
	
{
	
	if (v[i].charAt(0)=='-')
	
		C[i] = -(v[i].substring(1).hashCode());	
		
	else C[i] = (v[i].hashCode());
	
	
	
}



//System.out.println("taille X= "+X.length);
//System.out.println("taille Y= "+Y.length);

//System.out.println();

//System.out.println("X: "+s2);
//System.out.println();
//System.out.println("Y: "+s3);

//System.out.println();

Align al=new Align(X, Y, s4, s5);

if (show_time) 
	time=System.currentTimeMillis();
	

al.createAlign();	

//System.out.println("cout de l'alignement avant voisin "+al.cout);

if (vois) al.post_traitement_voisin_perte(C);

//System.out.println("cout de l'alignement apres voisin "+al.cout);

//al.Affichage();

if (dp) al.post_traitement_dupVperte();

if (ps) al.post_traitement_pertes_succ();

if (rp) al.pertes_succ();

if (show_time)

{
	time=System.currentTimeMillis()-time;

	String format = String.format("%%0%dd", 2);  
	  
	String seconds = String.format(format, (time/1000) % 60);  
	
	String minutes = String.format(format, ((time/1000) % 3600) / 60);  
	
	String hours = String.format(format, (time/1000) / 3600); 
	
	hhmmss=hours+":"+minutes+":"+seconds;

}



//al.get_distribution_events();  // determiner la distribution des �v�nements

cost=al.cout;

//dist_dup=al.stat_dup;

//dist_loss=al.stat_loss;

//dist_inver=al.stat_inver;

//dist_trans=al.stat_trans;

//System.out.println("cout de l'alignement apres pertes_succ "+al.cout);

//------------------- debut de changement

if (show_align) al.AlignOut1();

al.Affichage2();

if (show_time) System.out.println("Running time "+hhmmss);


//al.Affichage();

//System.out.println();

//System.out.println("Genome X = "+s2);

//System.out.println();

//System.out.println("Genome Y = "+s3);

//------------------- fin de changement


//debut de la construction de l'ancetre

createAncestor(al,vois);


}	




public static void createAncestor(Align al, boolean vois)

{

int Xlength=al.isCoveredX.length, Ylength=al.isCoveredY.length;

int i=0,j=0;

while (i<Xlength && j<Ylength)
	
	
	{
		if (al.isCoveredX[i].dup.value<=3)  // match, substitution, transposition ou inversion
	       
			{
			
			   if (al.isCoveredX[i].dup.value==2) // transposition
			
			   		{
				   
				       for (int h=al.isCoveredX[i].dup.y; h<=al.isCoveredX[i].dup.x; h++)
				    	   
				   		ancetre=ancetre+s4[h]+",";	
				   		
				       i=al.isCoveredX[i].dup.x+1;
				       
			   		}
			   
			   
				   
			   else   // match , substitution ou inversion
				   
				   
			   	{
				   
				   while (al.isCoveredY[j].dup.value!=0 && al.isCoveredY[j].dup.value!=1 && al.isCoveredY[j].dup.value!=3)
				   
				   		{
					   
					   //System.out.println("al.isCoveredY[j].dup.value = "+al.isCoveredY[j].dup.value);
					   
					   		if (al.isCoveredY[j].dup.value!=9)	  // Si ce n'est pas une perte dans X
					   			j=al.isCoveredY[j].dup.x+1; // ingnorer
					   			
					   		else
					   			
					   			{
					   			
					   			for (int h=al.isCoveredY[j].dup.y; h<=al.isCoveredY[j].dup.x; h++)
					   			
					   				ancetre=ancetre+s5[h]+",";
					   			
					   			
					   			j=al.isCoveredY[j].dup.x+1;
					   			
					   			}
					   		
					   		
					   			
				   		}
					
				   
				   // commencer le match la substitution ou l'inversion entre X et Y
				   
				   
				   for (int h=al.isCoveredX[i].dup.y; h<=al.isCoveredX[i].dup.x; h++)
					   
					   ancetre=ancetre+s4[h]+",";  // on l'a choisi dans X
				   
				   		i=al.isCoveredX[i].dup.x+1;
				   		j=al.isCoveredY[j].dup.x+1;
				   		
				   
			   	}
			
			   
			   
			
			}
		
		else
			
			{
			
			// commencer avec X
			
			 if (al.isCoveredX[i].dup.value==4)    
		
			 {
				 
				 for (int h=al.isCoveredX[i].dup.y; h<=al.isCoveredX[i].dup.x; h++)
					   
					   ancetre=ancetre+s4[h]+",";  // on l'a choisi dans X
			
			 }	 
				 
				 
			 i=al.isCoveredX[i].dup.x+1;
			
			 
			}
		
		
	}


//if (i<Xlength) System.out.println("Erreur");


while (i<Xlength)
	
	
{
	
	
		// commencer avec X
		
		 if (al.isCoveredX[i].dup.value==4)    
	
		 {
			 
			 for (int h=al.isCoveredX[i].dup.y; h<=al.isCoveredX[i].dup.x; h++)
				   
				   ancetre=ancetre+s4[h]+",";  // on l'a choisi dans X
		
		 }	 
			 
			 
		 i=al.isCoveredX[i].dup.x+1;
	
	
}


while (j<Ylength)

	{
	
	if (al.isCoveredY[j].dup.value==9)    
		
	 {
		 
		 for (int h=al.isCoveredY[j].dup.y; h<=al.isCoveredY[j].dup.x; h++)
			   
			   ancetre=ancetre+s5[h]+",";  // on l'a choisi dans X
	
	 }	 
		 
	 j=al.isCoveredY[j].dup.x+1;
	}
	
	
//System.out.println();	
System.out.println();


System.out.println(">Ancestor:");

ancetre=ancetre.substring(0, ancetre.length()-1);

System.out.println(ancetre);

System.out.println();	
System.out.println();


//System.out.println(">Distance : ");


//fin de la construction de l'ancetre



//System.out.println();
//System.out.println();


//if (vois) System.out.println("Post-traitement du voisinage : Oui");
//else System.out.println("Post-traitement du voisinage : Non");

}

	
}
	
	
	

	


	
