import QtQuick 2.12
import QtQuick.Controls 2.12
import QtQuick.Window 2.12

import com.example.MandelImageCombiner 1.0;
import com.example.MandelModel 1.0;
import com.example.JuliaModel 1.0;
import com.example.LaguerreModel 1.0;

Window {
    width: 640
    height: 480
    visible: true
    title: qsTr("Mandelbrot")
    //Keys.forwardTo: mouseArea //can't attach Keys to Window
    Component.onCompleted: { mandelModel.startRunning(); juliaModel.startRunning(); laguerreModel.startRunning(); }

    Column {
        anchors.left: parent.left;
        anchors.top: parent.top;
        id: toprow;
        //width: 150
        Text {
            id: labelXY;
            text: "fluff";
        }
        Text {
            id: labelInfoGen;
            text: "bluff";
        }
        Text {
            id: labelInfoSpec;
            text: "truff";
        }
    }

    Column {
        anchors.left: parent.left;
        anchors.right: rgroupView.left;
        anchors.top: parent.top;
        id: sliderCol;
        //width: 150
        Slider {
            anchors.right: parent.right;
            id: extAngleZoom;
            from: 0
            to: 10
            value: 0
            width: 50
            onValueChanged: { mandelModel.extAngleZoom=extAngleZoom.value; juliaModel.extAngleZoom=extAngleZoom.value; }
            ToolTip {
                text: "Angle zoom: "+extAngleZoom.value;
            }
        }
        CheckBox {
            anchors.right: parent.right;
            id: checkMandelPause
            checked: false
            text: "Pause"
            height: 13
            visible: rbuttonViewMand.checked
            onCheckedChanged: mandelModel.pause(checked);
        }
        CheckBox {
            anchors.right: parent.right;
            id: checkJuliaPause
            checked: false
            text: "Pause"
            height: 13
            visible: rbuttonViewJulia.checked
            onCheckedChanged: juliaModel.pause(checked);
        }
        CheckBox {
            anchors.right: parent.right;
            id: checkLaguPause
            checked: false
            text: "Pause"
            height: 13
            visible: rbuttonViewLagu.checked
            onCheckedChanged: laguerreModel.pause(checked);
        }
    }

    ButtonGroup { id: bgroupView }
    Column {
        anchors.right: comboColumn.left
        anchors.top: parent.top;
        id: rgroupView
        RadioButton {
            id: rbuttonViewMand
            checked: true
            text: "Mand"
            height: 13
            //implicitIndicatorHeight: 13
            ButtonGroup.group: bgroupView
        }
        RadioButton {
            id: rbuttonViewJulia
            checked: true
            text: "Julia"
            height: 13
            //implicitIndicatorHeight: 13
            ButtonGroup.group: bgroupView
        }
        RadioButton {
            id: rbuttonViewLagu
            text: "Lagu"
            height: 13
            //implicitIndicatorWidth: 13
            ButtonGroup.group: bgroupView
        }
    }

    Column {
        id: comboColumn
        anchors.top: parent.top;
        anchors.right: parent.right
        Row {
            anchors.right: parent.right
            Rectangle {
                id: busyIndicator
                width: 20
                height: 20
                radius: 10
                color: "blue"
            }
            ComboBox {
                textRole: "text"
                valueRole: "key"
                currentIndex: 1 //cls
                visible: rbuttonViewMand.checked
                model: ListModel {
                    id: mandel_paintstyle_model
                    ListElement { text: "Kind"; key: 0 }//mandelModel.paintStyleKind }
                    ListElement { text: "Classic"; key: 1 }//mandelModel.paintStyleCls }
                    ListElement { text: "Exter"; key: 2 }//mandelModel.paintStyleExter  }
                    ListElement { text: "ExtAng"; key: 3 }//mandelModel.paintStyleExter  }
                    ListElement { text: "Inter"; key: 4 }//mandelModel.paintStyleInter  }
                    ListElement { text: "Near"; key: 5 }//mandelModel.paintStyleNear  }
                    ListElement { text: "d/dz"; key: 6 }//mandelModel.paintStyleFZ  }
                    ListElement { text: "d/dc"; key: 7 }//mandelModel.paintStyleFC  }
                }
                onActivated: mandelModel.selectedPaintStyle=mandel_paintstyle_model.get(currentIndex).key;
            }
            ComboBox {
                textRole: "text"
                valueRole: "key"
                currentIndex: 1 //cls
                visible: rbuttonViewJulia.checked
                model: ListModel {
                    id: julia_paintstyle_model
                    ListElement { text: "Kind"; key: 0 }//mandelModel.paintStyleKind }
                    ListElement { text: "Classic"; key: 1 }//mandelModel.paintStyleCls }
                    ListElement { text: "Exter"; key: 2 }//mandelModel.paintStyleExter  }
                    ListElement { text: "ExtAng"; key: 3 }//mandelModel.paintStyleExter  }
                    ListElement { text: "Inter"; key: 4 }//mandelModel.paintStyleInter  }
                    ListElement { text: "Near"; key: 5 }//mandelModel.paintStyleNear  }
                    ListElement { text: "phi"; key: 6 }//
                    ListElement { text: "phi\'"; key: 7 }//
                    ListElement { text: "phi\'\'"; key: 8 }//
                    ListElement { text: "f\'"; key: 9 }//
                    ListElement { text: "|f'|<1"; key: 10 }//
                }
                onActivated: juliaModel.selectedPaintStyle=julia_paintstyle_model.get(currentIndex).key;
            }
            ComboBox {
                textRole: "text"
                valueRole: "key"
                currentIndex: 0 //cls
                visible: rbuttonViewLagu.checked
                model: ListModel {
                    id: lagu_paintstyle_model
                    ListElement { text: "Classic"; key: 0 }//laguModel.paintStyleCls
                    ListElement { text: "n-th f'"; key: 1 }//laguModel.paintStyleNthFz
                    ListElement { text: "n-th f'-1"; key: 2 }//laguModel.paintStyleNthFz1
                    ListElement { text: "roots"; key: 3 }//laguModel.paintStyleRoots
                }
                onActivated: laguerreModel.selectedPaintStyle=lagu_paintstyle_model.get(currentIndex).key;
            }
        }
        Row {
            anchors.right: parent.right
            Text {
                id: busyIndicator2
                text: "0"
            }
            ComboBox {
                textRole: "text"
                valueRole: "key"
                currentIndex: 0 //double
                model: ListModel {
                    id: precision_model
                    ListElement { text: "Double"; key: 0 }
                    ListElement { text: "Float128"; key: 1 }
                    ListElement { text: "DDouble"; key: 2 }
                    ListElement { text: "QDouble"; key: 3 }
                    ListElement { text: "Real642"; key: 4 }
                }
                onActivated:
                {
                    var prec=precision_model.get(currentIndex).key;
                    mandelModel.selectedPrecision=prec;
                    juliaModel.selectedPrecision=prec;
                    laguerreModel.selectedPrecision=prec;
                }
            }
        }
    }

    MouseArea {
        id: mouseArea
        anchors.left: parent.left;
        anchors.right: parent.right;
        anchors.top: toprow.bottom;
        anchors.bottom: parent.bottom;
        hoverEnabled: true;
        acceptedButtons: Qt.LeftButton | Qt.RightButton
        property int dragging: 0;
        property int drag_last_x;
        property int drag_last_y;
        Keys.onPressed: { //never triggers
            console.log(event);
            if (event.key===Qt.Key_E)
                mandelModel.paintOrbit(mandelImageCombiner.getOverlayImage(), mouse.x, mouse.y);
        }
        onPositionChanged:
        {
            if (dragging==1)
            {
                mandelModel.drag(mouse.x-drag_last_x, mouse.y-drag_last_y);
                drag_last_x=mouse.x;
                drag_last_y=mouse.y;
            };
            if (dragging==2)
            {
                juliaModel.drag(mouse.x-drag_last_x, mouse.y-drag_last_y);
                drag_last_x=mouse.x;
                drag_last_y=mouse.y;
            };
            if (dragging==3)
            {
                laguerreModel.drag(mouse.x-drag_last_x, mouse.y-drag_last_y);
                drag_last_x=mouse.x;
                drag_last_y=mouse.y;
            };
            //TODO: else copy from    ((mandel.mousePt.c.re<>view.orbit.re) or (mandel.mousePt.c.im<>view.orbit.im)) then
            //labelCX.text=mandelModel.pixelXtoRE_str(mouse.x);
            //labelCY.text=mandelModel.pixelYtoIM_str(mouse.y);
            if (rbuttonViewMand.checked)
              mandelModel.paintOrbit(mandelImageCombiner.getOverlayImage(), mouse.x, mouse.y);
            if (rbuttonViewJulia.checked)
              juliaModel.paintOrbit(juliaImageCombiner.getOverlayImage(), mouse.x, mouse.y);
            if (rbuttonViewLagu.checked)
            {
              laguerreModel.paintOrbit(laguerreImageCombiner.getOverlayImage(), mouse.x, mouse.y);
              if (laguerreModel==null)
                console.log("laguModel missing in action");
            }
            //labelCX.text=mouse.x;
            //labelCY.text=mouse.y;
            laguerreModel.getTimes();
        }
        onPressed: {
            if (mouse.button==Qt.LeftButton)
            {
                if (rbuttonViewMand.checked)
                  dragging=1;
                if (rbuttonViewJulia.checked)
                  dragging=2;
                if (rbuttonViewLagu.checked)
                  dragging=3;
                drag_last_x=mouse.x;
                drag_last_y=mouse.y;
            };
        }
        onReleased: {
            dragging=0;
            if (mouse.button==Qt.LeftButton)
            {
                if (rbuttonViewMand.checked)
                  mandelModel.paintOrbit(mandelImageCombiner.getOverlayImage(), mouse.x, mouse.y);
                if (rbuttonViewJulia.checked)
                  juliaModel.paintOrbit(mandelImageCombiner.getOverlayImage(), mouse.x, mouse.y);
                if (rbuttonViewLagu.checked)
                  laguerreModel.paintOrbit(laguerreImageCombiner.getOverlayImage(), mouse.x, mouse.y);
            }
            if (mouse.button==Qt.RightButton)
            {
                mainPopupMenu.popup(mouse.x, mouse.y);
            }
        }

        onWheel: {
            if (wheel.angleDelta.y>0)
            {
              if (rbuttonViewMand.checked)
                mandelModel.zoom(wheel.x, wheel.y, +1);
              else if (rbuttonViewJulia.checked)
                juliaModel.zoom(wheel.x, wheel.y, +1);
              else if (rbuttonViewLagu.checked)
                laguerreModel.zoom(wheel.x, wheel.y, +1);
            }
            else if (wheel.angleDelta.y<0)
            {
              if (rbuttonViewMand.checked)
                mandelModel.zoom(wheel.x, wheel.y, -1);
              else if (rbuttonViewJulia.checked)
                juliaModel.zoom(wheel.x, wheel.y, -1);
              else if (rbuttonViewLagu.checked)
                laguerreModel.zoom(wheel.x, wheel.y, -1);
            }
            wheel.accepted=true;
        }
        onHeightChanged: {
            mandelImageCombiner.resetBgImage(width, height, 0xff000000);
            mandelModel.setImageSize(width, height);
            mandelModel.writeToImage(mandelImageCombiner.getBaseImage());
            mandelImageCombiner.update();
            juliaImageCombiner.resetBgImage(width, height, 0xff000000);
            juliaModel.setImageSize(width, height);
            juliaModel.writeToImage(juliaImageCombiner.getBaseImage());
            juliaImageCombiner.update();
            laguerreImageCombiner.resetBgImage(width, height, 0xff000000);
            laguerreModel.setImageSize(width, height);
            laguerreModel.writeToImage(laguerreImageCombiner.getBaseImage());
            laguerreImageCombiner.update();
        }
        onWidthChanged: {
            mandelImageCombiner.resetBgImage(width, height, 0xff000000);
            mandelModel.setImageSize(width, height);
            mandelModel.writeToImage(mandelImageCombiner.getBaseImage());
            mandelImageCombiner.update();
            juliaImageCombiner.resetBgImage(width, height, 0xff000000);
            juliaModel.setImageSize(width, height);
            juliaModel.writeToImage(juliaImageCombiner.getBaseImage());
            juliaImageCombiner.update();
            laguerreImageCombiner.resetBgImage(width, height, 0xff000000);
            laguerreModel.setImageSize(width, height);
            laguerreModel.writeToImage(laguerreImageCombiner.getBaseImage());
            laguerreImageCombiner.update();
        }

        MandelImageCombiner {
            id: mandelImageCombiner;
            visible: rbuttonViewMand.checked
            Component.onCompleted: {
                resetBgImage(width, height, 0xff000000);
                mandelModel.setImageSize(width, height);
                //setBaseImage(mandelModel.getAsImage());
            }
            anchors.fill: parent;
        }

        MandelImageCombiner {
            id: juliaImageCombiner;
            visible: rbuttonViewJulia.checked
            Component.onCompleted: {
                resetBgImage(width, height, 0xff000000);
                juliaModel.setImageSize(width, height);
                //setBaseImage(mandelModel.getAsImage());
            }
            anchors.fill: parent;
        }

        MandelImageCombiner {
            id: laguerreImageCombiner;
            visible: rbuttonViewLagu.checked
            Component.onCompleted: {
                resetBgImage(width, height, 0xff000000);
                laguerreModel.setImageSize(width, height);
                //setBaseImage(mandelModel.getAsImage());
            }
            anchors.fill: parent;
        }

        Menu {
            id: mainPopupMenu
            /*MenuItem {
                text: "test1"
                onTriggered: console.log("clicked test1");
            }*/
            MenuItem {
                text: "set Lagu"
                onTriggered: { laguerreModel.setParams(mandelModel.viewInfo); juliaModel.setParams(mandelModel.viewInfo); }
            }
            Menu {
                title: "Presets (M)..."
                visible:rbuttonViewMand.checked
                Repeater {
                    model: [
                        {viewRe: -0.5, viewIm: 0, viewZoom: 1/128, caption: 'Main view'},
                        {viewRe: -0.75, viewIm: 0, viewZoom: 1/134217728, caption: 'Zoom on the 1-2 joint'},
                        {viewRe: -0.125, viewIm: 0.649519052836944556, viewZoom: 1/268435456, caption: 'Zoom on the 1-3 joint'},
                        {viewRe: -1.401155524, viewIm: 0, viewZoom: 1/8589934592, caption: 'A 1536 atom on real axis'},
                        {viewRe: -0.7458232657523781425, viewIm: 0.1659563046792491885, viewZoom: 1/1152921504606846976.0, caption: 'An extremely small 980 atom'},
                        {viewRe: -0.16400146484, viewIm: 1.0978469849, viewZoom: 1/262144, caption: 'Detail of ray 2/9'},
                        {viewRe: -0.54418945312, viewIm: 0.65454101562, viewZoom: 1/2048, caption: 'Detail of head of 2/5 bulb'},//  julia.c.re:=-0.53540039062,  julia.c.im:=-0.66320800781,//}
                        {viewRe: -0.706765832036462260364, viewIm: 0.3540325897413834077558, viewZoom: 1/1125899906842624, caption: 'Just some zoom of a 408 atom'},
                        {viewRe: -1.262303108533615159097935, viewIm: 0.383593549747416773254917, viewZoom: 1/70368744177664, caption: 'just a zoom'},
                        {viewRe: -1.67440967369022279335495, viewIm: 4.716562790887584498e-5, viewZoom: 1/137438953472, caption: 'Embedded Julia set'},
                        {viewRe: -0.7498572, viewIm: -0.0150818, viewZoom: 1/8388608, caption: 'Deep pools'},
                        {viewRe: 0.2860166666955662541362, viewIm: 0.0115374014484412084383, viewZoom: 1/2305843009213693952, caption: 'Spiral'},
                        {viewRe: 0.28601560169483100000, viewIm: 0.01153748599423500000, viewZoom: 1/17592186044416, caption: 'Dense 4-spiral'},
                        {viewRe: -1.39473140503329, viewIm: 0.00450421402574, viewZoom: 1/17179869184, caption: 'Star'},
                        {viewRe: -1.39473141148726511, viewIm: 0.00450422164786152379, viewZoom: 1/1125899906842624, caption: 'Circled mandel'},
                        {viewRe: -1.94098287909277016167, viewIm: 0.00064812392023092638912, viewZoom: 1/2199023255552, caption: 'Minijulia'},
                        {viewRe: -0.74364388703715904315, viewIm: 0.13182590420531251074, viewZoom: 1/8796093022208, caption: 'Wikipedia'},
                        //where did this come from {viewRe: -0.7499999970459612, viewIm: 0.0000514981111337, viewZoom: 1/549755813888, caption: 'Bulb 30502/61005'},
                        {viewRe: -0.7499999985575414, viewIm: 0.0000514972471137, viewZoom: 1/549755813888, caption: 'Bulb 30502/61005'},
                        {viewRe: -1498800068/8589934592.0, viewIm: 5685911464/8589934592.0, viewZoom: 1/8589934592.0, caption: 'Trouble at 11016'},
                    ]
                    MenuItem {
                        text: modelData.caption
                        onTriggered: mandelModel.setView_double(modelData.viewRe, modelData.viewIm, modelData.viewZoom);
                    }
                }
            }
            Menu {
                title: "Presets (J)..."
                visible:rbuttonViewJulia.checked
                Repeater {
                    model: [
                        {viewRe: -0.5, viewIm: 0, viewZoom: 1/128, cRe: 0, cIm: 0, period: 1, caption: 'Main view'},
                        {viewRe: -0.7342796903103590, viewIm: 0.1835422439035028, viewZoom: 1/4294967296.0, cRe: -0.7342796875163913, cIm: 0.1835422653239220, period: 1156, caption: 'Noisy zoom'},
                    ]
                    MenuItem {
                        text: modelData.caption
                        onTriggered: juliaModel.setParams(mandelModel.makeViewInfo(modelData));
                    }
                }
            }
        }
    }



    MandelModel {
        id: mandelModel
    }

    JuliaModel {
        id: juliaModel
    }

    LaguerreModel {
        id: laguerreModel
    }

    Timer {
        interval: 100;
        repeat: true;
        running: true;
        onTriggered: {
            if (rbuttonViewMand.checked)
            {
              mandelModel.writeToImage(mandelImageCombiner.getBaseImage());
              mandelImageCombiner.update();
              labelXY.text=mandelModel.getTextXY();
              labelInfoGen.text=mandelModel.getTextInfoGen();
              labelInfoSpec.text=mandelModel.getTextInfoSpec();
              mandelModel.getTimes();
              if (mandelModel.threadsWorking===0)
                busyIndicator.color="green";
              else if (mandelModel.threadsWorking<mandelModel.threadsMax)
                busyIndicator.color="yellow";
              //else if (mandelModel.threadsWorking>mandelModel.threadsMax)
              //  busyIndicator.color="blue"; //up to 2*threads during dragging... or more
              else
                busyIndicator.color="red";
              busyIndicator2.text=mandelModel.threadsWorking
            };
            if (rbuttonViewJulia.checked)
            {
              juliaModel.writeToImage(juliaImageCombiner.getBaseImage());
              juliaImageCombiner.update();
              labelXY.text=juliaModel.getTextXY();
              labelInfoGen.text=juliaModel.getTextInfoGen();
              labelInfoSpec.text=juliaModel.getTextInfoSpec();
              juliaModel.getTimes();
              if (juliaModel.threadsWorking===0)
                busyIndicator.color="green";
              else if (juliaModel.threadsWorking<juliaModel.threadsMax)
                busyIndicator.color="yellow";
              //else if (mandelModel.threadsWorking>mandelModel.threadsMax)
              //  busyIndicator.color="blue"; //up to 2*threads during dragging... or more
              else
                busyIndicator.color="red";
              busyIndicator2.text=juliaModel.threadsWorking
            };
            if (rbuttonViewLagu.checked)
            {
              laguerreModel.getTimes();
              laguerreModel.writeToImage(laguerreImageCombiner.getBaseImage());
              laguerreImageCombiner.update();
              labelXY.text=laguerreModel.getTextXY();
              labelInfoGen.text=laguerreModel.getTextInfoGen();
              labelInfoSpec.text=laguerreModel.getTextInfoSpec();
              busyIndicator.color="black";
              busyIndicator2.text="";
            };
        }
    }
}
