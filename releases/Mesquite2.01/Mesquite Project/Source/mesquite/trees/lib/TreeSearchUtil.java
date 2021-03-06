/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.
Version 2.01, December 2007.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
*/
package mesquite.trees.lib;

import mesquite.lib.*;
import mesquite.lib.duties.*;

public class TreeSearchUtil {

	/*.................................................................................................................*/
	private static void cleanUpSearch(MesquiteTimer timer, MesquiteNumber currentScore, MesquiteTree swapTree, AdjustableTree tree, int node, MesquiteModule ownerModule, MesquiteLong count,  boolean liveUpdates, boolean notify) {
		double seconds = timer.timeSinceLastInSeconds();
		MesquiteMessage.println("\n\nSearch completed.");
		MesquiteMessage.println("    Final score: " + currentScore.toString());
		MesquiteMessage.println("    Total number of rearrangements: " + count.toString());
		MesquiteMessage.println("    Time taken: " + seconds + " seconds");
		swapTree.standardize(node,true, false);
		tree.setToClone(swapTree);
		if (notify && tree instanceof Listened && !liveUpdates) {
			((Listened)tree).notifyListeners(ownerModule, new Notification(MesquiteListener.BRANCHES_REARRANGED));
		}
	}
	
	
	/*.................................................................................................................*/
	public static boolean tryRearrangement(long i, long total,MesquiteTimer timer, MesquiteNumber currentScore, MesquiteModule ownerModule, MesquiteTree swapTree, MesquiteTree tempTree, AdjustableTree tree, int node, TreeSwapper swapTask, NumberForTree numberTask, ProgressIndicator progIndicator, MesquiteLong count, MesquiteBoolean foundBetter, boolean smallerIsBetter, boolean liveUpdates, boolean notify) {
		count.increment();
		MesquiteNumber tempScore = new MesquiteNumber();
		if (progIndicator != null && i % 20==0) {
			if (progIndicator.isAborted()) {
					progIndicator.goAway();
					cleanUpSearch(timer, currentScore, swapTree,tree,node, ownerModule, count, liveUpdates, notify);
					return false;
			}
			progIndicator.spin();
		}
		
		MesquiteString rs =new MesquiteString("");
		swapTask.rearrange(tempTree, node, i);  //do the rearrangement!
		numberTask.calculateNumber(tempTree, tempScore, rs);
		if (i%100==0) {
			MesquiteMessage.print(".");
		}
		if (total%20==0) {
			CommandRecord.tick(numberTask.getName()+ ": " + currentScore.toString() + "   ("+total+" rearrangements)");
			if (progIndicator!=null)
				progIndicator.setText(numberTask.getName()+ ": " + currentScore.toString() + "   ("+total+" rearrangements)");
		}
		
		if ((smallerIsBetter && tempScore.isLessThan(currentScore)) || (!smallerIsBetter && tempScore.isMoreThan(currentScore))){
			currentScore.setValue(tempScore);
			swapTree.setToClone(tempTree);
			foundBetter.setValue(true);
			if (notify && tree instanceof Listened && liveUpdates) {
				swapTree.standardize(node,true, false);
				tree.setToClone(swapTree);
				((Listened)tree).notifyListeners(ownerModule, new Notification(MesquiteListener.BRANCHES_REARRANGED));
			}
			return true;
		} else
			tempTree.setToClone(swapTree);			
		return true;
	}
	
	/*.................................................................................................................*/
	public  static boolean searchForBetterTree(MesquiteModule ownerModule, AdjustableTree tree, int node, TreeSwapper swapTask, NumberForTree numberTask, RandomBetween rng, MesquiteString resultString, boolean smallerIsBetter, boolean liveUpdates, boolean notify){
		numberTask.initialize(tree);
		MesquiteTimer timer = new MesquiteTimer();
		timer.start();
		
		MesquiteString rs =new MesquiteString();
		MesquiteNumber currentScore = new MesquiteNumber();
		rs.setValue("");
		MesquiteTree swapTree = new MesquiteTree(tree.getTaxa());
		MesquiteTree tempTree = new MesquiteTree(tree.getTaxa());
		swapTree = tree.cloneTree();
		numberTask.calculateNumber(swapTree, currentScore, rs);
		MesquiteBoolean foundBetter= new MesquiteBoolean(true);
		MesquiteDialog.hideWizardForCalculation();
		
		ProgressIndicator progIndicator=null;
		if (ownerModule!=null && !MesquiteThread.isScripting())
			progIndicator=new ProgressIndicator(ownerModule.getProject(),ownerModule.getName(), "Searching for a better tree", 0, true);
		if (progIndicator!=null){
			progIndicator.setButtonMode(ProgressIndicator.OFFER_CONTINUE);
			progIndicator.setOfferContinueMessageString("Are you sure you want to stop the search?");
			progIndicator.start();
		}

		MesquiteLong count = new MesquiteLong(0);
		long total = 0;

		
		while(foundBetter.getValue()) {  // loop for improving tree
			tempTree.setToClone(swapTree);
			foundBetter.setValue(false);
			long numRearrangements = swapTask.numberOfRearrangements(swapTree, node);
			MesquiteMessage.print("\n  " + numberTask.getName()+ ": " + currentScore.toString()+ " ");
			
			int motherNode = tree.motherOfNode(node); //store information so we can find our original node even if it has been renumbered.
			int aTerminalOfNode=tree.leftmostTerminalOfNode(node);

			/* Rather than going through the rearrangements in order (which will mean that one generally starts in the same part of the tree, and thus one will cover that
			 * territory extremely well, to the detriment of other areas of the tree), we here pick a random rearrangement to serve as our starting point.  This is called "boundary".
			 * We then either go down from there to rearrangement zero, followed by going up from boundary to numRearrangements, or we do the reverse, first going up from boundary.
			 * This insures that we start in a random part of the tree each time.  This was a method suggested many years ago by DRM to DLS. The other method (going through rearrangements in order) 
			 * is deterministic; this method with a random boundary, however, will give a different search path each time it is run.  
			 * I did two comparisons of the two approaches, one with 27 taxa, the other with 53 taxa.
			 * 27 taxa:  Standard order of rearrangement:  80.7, 80.6, and 80.5 seconds.
			 *                with random starting point:   59.2, 27.4, 36.6, and 48.1 seconds.     
			 * 53 taxa:  Standard order of rearrangement:  1505 seconds.
			 *                with random starting point:   322 seconds.     
			 * Ideally one might want to use a more
			 * intelligent way to pick the boundary, but as that depends upon the nature of the objective function, and this code knows nothing but its value, that would be very hard.  */

			/* if one used the standard order, this is what the next section would look like:
			
			 for (long i=0; i<numRearrangements; i++) {
				if (!tryRearrangement( i,  timer, currentScore, ownerModule,  swapTree,  tempTree,  tree,  node, swapTask,  numberTask,  progIndicator,  count,  foundBetter,  smallerIsBetter,  liveUpdates,  notify))
					return false;
				if (foundBetter.getValue()) break;
			}
			*/

			CommandRecord.tick(numberTask.getName()+ ": " + currentScore.toString());
			if (progIndicator!=null)
				progIndicator.setText(numberTask.getName()+ ": " + currentScore.toString());
			long boundary=rng.randomLongBetween(0, numRearrangements-1);
			boolean downFirst  = rng.randomIntBetween(0, 1)>0;
			
			if (downFirst) {
				for (long i=boundary; i<numRearrangements; i++) {
					total++;
					if (!tryRearrangement( i,  total, timer, currentScore, ownerModule,  swapTree,  tempTree,  tree,  node, swapTask,  numberTask,  progIndicator,  count,  foundBetter,  smallerIsBetter,  liveUpdates,  notify))
						return false;
					if (foundBetter.getValue()) break;
				}
				if (!foundBetter.getValue()){  // then let's try the other rearrangements
					for (long i=boundary-1; i>=0; i--) {
						total++;
						if (!tryRearrangement( i,  total, timer, currentScore, ownerModule,  swapTree,  tempTree,  tree,  node, swapTask,  numberTask,  progIndicator,  count,  foundBetter,  smallerIsBetter,  liveUpdates,  notify))
							return false;
						if (foundBetter.getValue()) break;
					}
				}
			}
			else {
				for (long i=boundary-1; i>=0; i--) {
					total++;
					if (!tryRearrangement( i,  total, timer, currentScore, ownerModule,  swapTree,  tempTree,  tree, node, swapTask,  numberTask,  progIndicator,  count,  foundBetter,  smallerIsBetter,  liveUpdates,  notify))
						return false;
					if (foundBetter.getValue()) break;
				}
				if (!foundBetter.getValue()){  // then let's try the other rearrangements
					total++;
					for (long i=boundary; i<numRearrangements; i++) {
						if (!tryRearrangement( i,  total, timer, currentScore, ownerModule,  swapTree,  tempTree,  tree,  node, swapTask,  numberTask,  progIndicator,  count,  foundBetter,  smallerIsBetter,  liveUpdates,  notify))
							return false;
						if (foundBetter.getValue()) break;
					}
				}
			}
	
			
			for (int d = tree.firstDaughterOfNode(motherNode); tree.nodeExists(d); d = tree.nextSisterOfNode(d)) {
				if (tree.descendantOf(aTerminalOfNode, d)) { // then this is our original node
					node = d;
					break;
				}

			}
			
		}  // end while loop
		
		cleanUpSearch(timer, currentScore, swapTree,tree, node, ownerModule, count, liveUpdates, notify);

		if (progIndicator!=null)
			progIndicator.goAway();

		return true;
	}


}
