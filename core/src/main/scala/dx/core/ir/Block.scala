package dx.core.ir

import dx.util.Enum

import scala.reflect.ClassTag

/**
  * The kind of block input variable being referenced. A Computed
  * input is one that is computed from other inputs - currently
  * this is only used for the scatter variable.
  */
object InputKind extends Enum {
  type InputKind = Value
  val Required, Computed, Optional = Value
}

/**
  * These are the kinds of blocks that are run by the workflow-fragment-runner.
  * A block can have expressions, input ports, and output ports in the beginning.
  * The block can be one of these types:
  *   - Purely expressions (no asynchronous calls at any nesting level).
  *   - Call
  *      - with no subexpressions to evaluate
  *      - with subexpressions requiring evaluation
  *      - Fragment with expressions and one call
  *   - Conditional block
  *      - with exactly one call
  *      - a complex subblock
  *   - Scatter block
  *      - with exactly one call
  *      - with a complex subblock
  */
object BlockKind extends Enum {
  type BlockKind = Value
  val ExpressionsOnly, CallDirect, CallWithSubexpressions, CallFragment, ConditionalOneCall,
      ConditionalComplex, ScatterOneCall, ScatterComplex = Value
}

trait Block[Self <: Block[Self]] { this: Self =>

  /**
    * The kind of block this is.
    */
  def kind: BlockKind.BlockKind

  /**
    * The index of this block within its parent list.
    */
  def index: Int

  /**
    * Create a human readable name for a block of statements
    * 1. Ignore all private variables
    * 2. If there is a scatter/if, use that
    * 3. if there is at least one call, use the first one.
    * 4. If the entire block is made up of expressions, return None
    * @return
    */
  def getName: Option[String]

  def getSubBlock(index: Int): Self

  def getSubBlockAt(path: Vector[Int])(implicit tag: ClassTag[Self]): Self = {
    path.foldLeft(this) {
      case (subBlock: Self, index) => subBlock.getSubBlock(index)
    }
  }

  def inputNames: Set[DxName]

  def outputNames: Set[DxName]

  def prettyFormat: String
}

object Block {
  def getSubBlockAt[B <: Block[B]: ClassTag](topBlocks: Vector[B], path: Vector[Int]): B = {
    topBlocks(path.head).getSubBlockAt(path.tail)
  }
}
