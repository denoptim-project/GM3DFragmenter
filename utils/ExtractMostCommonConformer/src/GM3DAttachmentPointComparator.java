import java.util.Comparator;

/**
 *
 * @author Marco Foscato (University of Bergen) Discussions: Giovanni Occhipinti and Vidar R. Jensen (University of Bergen)
 */

public class GM3DAttachmentPointComparator implements Comparator<GM3DAttachmentPoint>
{
   @Override
   public int compare(GM3DAttachmentPoint ap1, GM3DAttachmentPoint ap2)
    {
        return ap1.compareTo(ap2);
    }

}

